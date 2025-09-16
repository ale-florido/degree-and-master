import os
from dotenv import load_dotenv

# Carga las variables de entorno desde un archivo .env
load_dotenv()

# Configuración de DB - Obtiene cada variable de entorno
DB_HOST = os.getenv("DB_HOST")        # Host de la base de datos
DB_PORT = os.getenv("DB_PORT")        # Puerto de la base de datos
DB_NAME = os.getenv("DB_NAME")        # Nombre de la base de datos
DB_USER = os.getenv("DB_USER")        # Usuario de la base de datos
DB_PASS = os.getenv("DB_PASS")        # Contraseña de la base de datos
DB_TYPE = os.getenv("DB_TYPE")        # Tipo de base de datos (mysql, postgres, etc.)
