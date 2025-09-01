const mysql = require('mysql2/promise');
// Importa el cliente MySQL compatible con promesas.

const HOST = 'mysql-sakila-test';
// Nombre del host del contenedor MySQL (debe estar en la misma red de Docker).

const USER = 'username';
// Usuario de la base de datos.

const PASSWORD = 'password';
// Contraseña del usuario de la base de datos.

const DATABASE = 'sakila';
// Nombre de la base de datos a consultar.

const PORT = 3306;
// Puerto por defecto de MySQL.

let sqlSchemaPrompt = '';
// Variable donde se almacenará el resultado final.

async function getSQLPrompt() {
    let connection;

    try {
        connection = await mysql.createConnection({
            host: HOST,
            user: USER,
            password: PASSWORD,
            database: DATABASE,
            port: PORT,
            connectTimeout: 5000 // el cliente MySQL esperará hasta 5 segundos para conectarse al servidor antes de fallar con un error de timeout
        });
        // Establece la conexión a la base de datos.

        // Consulta para obtener todas las tablas de la base de datos.
        const [tablesResult] = await connection.query(`
            SELECT table_name
            FROM information_schema.tables
            WHERE table_schema = ? AND table_type = 'BASE TABLE'
            ORDER BY table_name ASC
        `, [DATABASE]);

        if (tablesResult.length === 0) {
            // Si no hay tablas, guarda un mensaje de advertencia.
            sqlSchemaPrompt = `⚠️ No se encontraron tablas en la base de datos \`${DATABASE}\`.`;
            return;
        }

        // Encabezado del resultado.
        sqlSchemaPrompt = `✅ **Lista de tablas y atributos en la base de datos \`${DATABASE}\`:**\n\n`;

        for (const tableRow of tablesResult) {
            const tableName = tableRow.table_name || tableRow.TABLE_NAME;
            // Obtiene el nombre de la tabla.

            sqlSchemaPrompt += `- \`${tableName}\`\n`;

            // Consulta para obtener las columnas y tipos de datos de la tabla actual.
            const [columnsResult] = await connection.query(`
                SELECT column_name, data_type
                FROM information_schema.columns
                WHERE table_schema = ? AND table_name = ?
                ORDER BY ordinal_position ASC
            `, [DATABASE, tableName]);

            if (columnsResult.length > 0) {
                sqlSchemaPrompt += `  - **Atributos:**\n`;
                for (const col of columnsResult) {
                    const columnName = col.column_name || col.COLUMN_NAME || 'undefined';
                    const dataTypeRaw = col.data_type || col.DATA_TYPE;
                    const dataType = dataTypeRaw ? dataTypeRaw.toUpperCase() : 'UNKNOWN';
                    sqlSchemaPrompt += `    - \`${columnName}\` (${dataType})\n`;
                }
            } else {
                sqlSchemaPrompt += `  - *(Sin atributos encontrados)*\n`;
            }
        }

    } catch (err) {
        // Si ocurre un error, lo muestra y lo guarda en la variable de resultado.
        console.error('❌ Error en getSQLPrompt:', err);
        sqlSchemaPrompt = `❌ Error generando listado de tablas: ${err.message}`;
    } finally {
        // Cierra la conexión si está abierta.
        if (connection) await connection.end();
    }
}

// Ejecución principal del script.
(async () => {
    await getSQLPrompt();
    // Imprime el resultado final en consola (stdout), que será capturado por server.js.
    console.log(sqlSchemaPrompt);
})();
