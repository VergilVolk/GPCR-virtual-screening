FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PACER_M4_ROOT=/app \
    PACER_M4_API_ALLOW_EXECUTION=0

WORKDIR /app

COPY pyproject.toml README.md LICENSE ./
COPY pacer_m4 ./pacer_m4
COPY project ./project

RUN python -m pip install --no-cache-dir ".[api,core]"

EXPOSE 8000

HEALTHCHECK --interval=30s --timeout=5s --start-period=10s --retries=3 \
  CMD python -c "import urllib.request; urllib.request.urlopen('http://127.0.0.1:8000/health', timeout=3)"

CMD ["uvicorn", "pacer_m4.api:app", "--host", "0.0.0.0", "--port", "8000"]

