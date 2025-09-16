from sqlalchemy import create_engine, text
from sqlalchemy.exc import SQLAlchemyError
from .config import DB_TYPE, DB_USER, DB_PASS, DB_HOST, DB_PORT, DB_NAME

# Construye la URL de conexión a la base de datos usando el formato de SQLAlchemy
# Formato: mysql+mysqlconnector://usuario:contraseña@host:puerto/nombre_base_datos
DATABASE_URL = f"mysql+mysqlconnector://{DB_USER}:{DB_PASS}@{DB_HOST}:{DB_PORT}/{DB_NAME}"

# Crea el motor de SQLAlchemy que maneja las conexiones con la base de datos
# El motor sirve como pool de conexiones y punto centralizado para operaciones DB
engine = create_engine(DATABASE_URL)

def execute_sql(query: str, params: dict = None):
    """Ejecuta consultas SQL de lectura solamente"""
    
    # Comandos permitidos
    READ_COMMANDS = {"select", "show", "describe", "explain"}
    first_word = query.strip().split()[0].lower()
    # Extrae la primera palabra del query (en minúsculas) para validación
    # Ejemplo: "SELECT * FROM users" -> "select"
    
    if first_word not in READ_COMMANDS:
        return {"error": "Solo consultas de lectura permitidas"}
    
    try:
        # Establece una conexión con la base de datos usando context manager
        with engine.connect() as connection:
            result = connection.execute(text(query), params or {})
            # Convierte cada fila del resultado en un diccionario
            return [dict(row) for row in result.mappings()]
    except SQLAlchemyError as e:
        return {"error": str(e)}