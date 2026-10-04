# ChronoCell-5D 4.0: the app, the test suite and the validation harness in one CPU image.
#   docker build -t chronocell .
#   docker run --rm -p 8501:8501 chronocell                              # the app on http://localhost:8501
#   docker run --rm chronocell python -m pytest -q                       # the test suite
#   docker run --rm -v "$PWD/out:/app/paper/figures" chronocell python -m validation.reproduce   # paper figures
FROM python:3.10-slim

ENV PYTHONDONTWRITEBYTECODE=1     PYTHONUNBUFFERED=1     PIP_NO_CACHE_DIR=1     PIP_DISABLE_PIP_VERSION_CHECK=1

WORKDIR /app
COPY requirements.lock requirements-validation.txt ./
RUN pip install --index-url https://download.pytorch.org/whl/cpu torch==2.11.0  && pip install -r requirements.lock -r requirements-validation.txt  && pip freeze > /app/pip-freeze.txt

COPY . .
RUN useradd --create-home chronocell && chown -R chronocell /app
USER chronocell

EXPOSE 8501
HEALTHCHECK --interval=30s --timeout=5s CMD python -c "import urllib.request; urllib.request.urlopen('http://localhost:8501/_stcore/health')"
CMD ["streamlit", "run", "app.py", "--server.address=0.0.0.0", "--server.port=8501", "--server.headless=true"]
