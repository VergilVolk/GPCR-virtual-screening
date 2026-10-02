FROM python:3.9.23-slim-bookworm

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PACER_M4_ROOT=/app \
    PACER_M4_API_ALLOW_EXECUTION=0

WORKDIR /app

COPY requirements.txt requirements-api.txt pyproject.toml README.md LICENSE ./
COPY pacer_m4 ./pacer_m4
COPY project ./project

# Install the CPU wheel first.  The exact core requirements then keep the
# already-satisfied torch build instead of pulling CUDA runtime packages.
RUN python -m pip install --no-cache-dir --upgrade pip \
    && python -m pip install --no-cache-dir \
       torch==2.8.0 --index-url https://download.pytorch.org/whl/cpu \
    && python -m pip install --no-cache-dir -r requirements.txt \
    && python -m pip install --no-cache-dir -r requirements-api.txt \
    && python -m pip install --no-cache-dir --no-deps . \
    && python project/scripts/audit_core_environment.py \
    && python -m pacer_m4 capabilities --output /tmp/pacer_m4_capabilities.json

EXPOSE 8000

HEALTHCHECK --interval=30s --timeout=5s --start-period=10s --retries=3 \
  CMD python -c "import urllib.request; urllib.request.urlopen('http://127.0.0.1:8000/health', timeout=3)"

CMD ["uvicorn", "pacer_m4.api:app", "--host", "0.0.0.0", "--port", "8000"]
