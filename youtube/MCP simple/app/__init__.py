from fastmcp import FastMCP
from .db_client import execute_sql

# Crea una instancia de FastMCP con el nombre "MCP SQL Server"
# Este nombre identifica el servidor MCP para los clientes que se conecten
mcp = FastMCP("MCP SQL Server")

@mcp.tool
def execute_sql_tool(query: str, params: dict = None):
    """Ejecuta una consulta SQL de lectura"""
    return execute_sql(query, params)
