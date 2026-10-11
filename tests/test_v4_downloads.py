"""Downloads of the predictor's inputs (chronocell.predict): the UCSC chromosome sequence and ENCODE CTCF
peaks. A local HTTP server stands in for UCSC and ENCODE (no internet in tests); it honours Range
requests, can cut a response short, and can serve a wrong MD5."""

from __future__ import annotations

import gzip
import hashlib
import json
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

import pytest

from chronocell import predict as PD

PAYLOAD = gzip.compress(b">chr21\n" + b"ACGTTGCA" * 200_000 + b"\n")       # ~ 3 kB compressed, 1.6 MB of sequence


class _Server:
    """Files by path; options: cut the next n responses after k bytes, or ignore Range headers."""

    def __init__(self, files: dict[str, bytes]):
        self.files = files
        self.requests: list[tuple[str, str | None]] = []
        self.cut_next = 0
        self.cut_at = 1000
        srv = self

        class H(BaseHTTPRequestHandler):
            def log_message(self, *a):            # quiet
                pass

            def do_GET(self):
                rng = self.headers.get("Range")
                srv.requests.append((self.path, rng))
                body = srv.files.get(self.path)
                if body is None:
                    self.send_response(404)
                    self.end_headers()
                    return
                start = 0
                if rng and rng.startswith("bytes="):
                    start = int(rng[6:].split("-")[0])
                    self.send_response(206)
                    self.send_header("Content-Range", f"bytes {start}-{len(body) - 1}/{len(body)}")
                else:
                    self.send_response(200)
                part = body[start:]
                self.send_header("Content-Length", str(len(part)))
                self.end_headers()
                if srv.cut_next > 0 and self.path.endswith(".gz"):   # promise the whole length, send less
                    srv.cut_next -= 1
                    self.wfile.write(part[: srv.cut_at])
                    self.wfile.flush()
                    self.close_connection = True
                    return
                self.wfile.write(part)

        self.httpd = ThreadingHTTPServer(("127.0.0.1", 0), H)
        self.url = f"http://127.0.0.1:{self.httpd.server_address[1]}"
        threading.Thread(target=self.httpd.serve_forever, daemon=True).start()

    def close(self):
        self.httpd.shutdown()
        self.httpd.server_close()


@pytest.fixture()
def server(monkeypatch):
    md5 = hashlib.md5(PAYLOAD).hexdigest()
    s = _Server({"/hg38/chromosomes/chr21.fa.gz": PAYLOAD,
                 "/hg38/chromosomes/md5sum.txt": f"{md5}  chr21.fa.gz\n0123456789abcdef0123456789abcdef  chr22.fa.gz\n".encode(),
                 "/files/ENCFFTEST/@@download/ENCFFTEST.bed.gz": PAYLOAD})
    monkeypatch.setattr(PD, "UCSC_CHROMOSOMES", s.url + "/{assembly}/chromosomes/")
    monkeypatch.setattr(PD, "ENCODE_DOWNLOAD", s.url + "/files/{acc}/@@download/{acc}.bed.gz")
    monkeypatch.setattr(time, "sleep", lambda _s: None)                          # retries without waiting
    yield s
    s.close()


def test_sequence_downloads_and_is_checked_against_ucsc_md5(server, tmp_path):
    seen = []
    out = PD.fetch_chromosome_fasta("hg38", "chr21", tmp_path, lambda d, t: seen.append((d, t)))
    assert out.read_bytes() == PAYLOAD and out.name == "hg38_chr21.fa.gz"
    assert seen[-1] == (len(PAYLOAD), len(PAYLOAD))
    assert PD.read_fasta(out) == b"ACGTTGCA" * 200_000
    n = len(server.requests)
    assert PD.fetch_chromosome_fasta("hg38", "chr21", tmp_path) == out and len(server.requests) == n   # cached


def test_interrupted_download_resumes_with_a_range_request(server, tmp_path):
    part = tmp_path / "hg38_chr21.fa.gz.part"
    part.write_bytes(PAYLOAD[:700])                                 # left over from an interrupted run
    out = PD.fetch_chromosome_fasta("hg38", "chr21", tmp_path)
    assert out.read_bytes() == PAYLOAD and not part.exists()
    assert ("/hg38/chromosomes/chr21.fa.gz", "bytes=700-") in server.requests


def test_short_read_is_retried_and_resumed(server, tmp_path):
    server.cut_next, server.cut_at = 2, 500                         # two responses end early
    out = PD.fetch_chromosome_fasta("hg38", "chr21", tmp_path)
    assert out.read_bytes() == PAYLOAD
    ranges = [r for p, r in server.requests if p.endswith("chr21.fa.gz")]
    assert ranges == [None, "bytes=500-", "bytes=1000-"]


def test_md5_mismatch_deletes_the_file(server, tmp_path):
    server.files["/hg38/chromosomes/md5sum.txt"] = b"0" * 32 + b"  chr21.fa.gz\n"
    with pytest.raises(OSError, match="MD5"):
        PD.fetch_chromosome_fasta("hg38", "chr21", tmp_path)
    assert list(tmp_path.iterdir()) == []                           # neither the .part nor the final file


def test_no_published_md5_means_no_download(server, tmp_path):
    with pytest.raises(OSError, match="unverifiable"):
        PD.fetch_chromosome_fasta("hg38", "chr7", tmp_path)
    assert not any(p.endswith("chr7.fa.gz") for p, _ in server.requests)


def test_encode_peaks_by_accession_with_md5_and_resume(server, tmp_path):
    md5 = hashlib.md5(PAYLOAD).hexdigest()
    (tmp_path / "ENCFFTEST.bed.gz.part").write_bytes(PAYLOAD[:300])
    server.cut_next = 1
    out = PD.fetch_encode_peaks("ENCFFTEST", md5, tmp_path)
    assert out.read_bytes() == PAYLOAD
    ranges = [r for p, r in server.requests if "ENCFFTEST" in p]
    assert ranges[0] == "bytes=300-" and len(ranges) == 2
    with pytest.raises(OSError, match="MD5"):
        PD.fetch_encode_peaks("ENCFFTEST", "f" * 32, tmp_path / "again")
    assert not (tmp_path / "again" / "ENCFFTEST.bed.gz").exists()
    assert not (tmp_path / "again" / "ENCFFTEST.bed.gz.part").exists()


def test_encode_sources_are_the_listed_gate5_inputs():
    srcs = {s["cell_line"]: s for s in PD.encode_ctcf_sources()}
    assert set(srcs) == {"IMR90", "A549", "K562", "HCT116"}
    rows = {r["key"]: r for r in json.loads(PD.SOURCES_PATH.read_text(encoding="utf-8"))}
    assert srcs["IMR90"]["accession"] == "ENCFF670ULH" and srcs["IMR90"]["md5"] == rows["encode_ctcf_imr90"]["md5"]
    assert all(len(s["md5"]) == 32 and s["assembly"] == "GRCh38" and s["experiment"].startswith("ENCSR")
               for s in srcs.values())
