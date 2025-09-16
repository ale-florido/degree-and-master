# Importa la instancia del servidor MCP desde el módulo app
# Esto asume que en app/__init__.py hay una variable 'mcp' exportada
from app import mcp

if __name__ == "__main__":
    mcp.run(
        transport="sse", # SSE es ideal para comunicación unidireccional del servidor al cliente
        host="0.0.0.0", # Escucha en todas las interfaces de red (0.0.0.0)
        port=8000, # Puerto donde el servidor estará escuchando
        log_level="info",  # Nivel de logging: "info" muestra información general de operación
        path="/mcp" # Ruta base para los endpoints del servidor MCP
    )