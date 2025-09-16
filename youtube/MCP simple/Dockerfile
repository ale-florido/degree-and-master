FROM python:3.11-slim

WORKDIR /app

# Instalar dependencias para MySQL y compilación
RUN apt-get update && \
    apt-get install -y default-libmysqlclient-dev build-essential && \
    rm -rf /var/lib/apt/lists/*

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY app/ ./app/
COPY .env ./

EXPOSE 8000

# Ejecutar el servidor MCP (FastMCP)
CMD ["python", "-m", "app.run_mcp"]




