const mysql = require('mysql2/promise');

// Configuración de conexión
const HOST = 'mysql-sakila-test';
const USER = 'root';
const PASSWORD = 'password';
const DATABASE = 'sakila'; //NOMBRE DE sakila SI ES PARA LA CUSTOM FUNCITON DE sakila, O world SI ES PARA LA CUSTOM FUNCTION DE world
const PORT = 3306;

let sqlSchemaPrompt = '';

async function getSQLPrompt() {
    let connection;

    try {
        connection = await mysql.createConnection({
            host: HOST,
            user: USER,
            password: PASSWORD,
            database: DATABASE,
            port: PORT,
            connectTimeout: 5000
        });

        // Obtener todos los nombres de tablas
        const [tablesResult] = await connection.query(`
            SELECT table_name
            FROM information_schema.tables
            WHERE table_schema = ? AND table_type = 'BASE TABLE'
            ORDER BY table_name ASC
        `, [DATABASE]);

        if (tablesResult.length === 0) {
            sqlSchemaPrompt = `⚠️ No se encontraron tablas en la base de datos \`${DATABASE}\`.`;
            return;
        }

        sqlSchemaPrompt = `✅ **Lista de tablas y atributos en la base de datos \`${DATABASE}\`:**\n\n`;

        for (const tableRow of tablesResult) {
            const tableName = tableRow.table_name || tableRow.TABLE_NAME;

            sqlSchemaPrompt += `- \`${tableName}\`\n`;

            // Obtener columnas y tipos de dato de la tabla actual
            const [columnsResult] = await connection.query(`
                SELECT column_name, data_type
                FROM information_schema.columns
                WHERE table_schema = ? AND table_name = ?
                ORDER BY ordinal_position ASC
            `, [DATABASE, tableName]);

            // DEBUG: ver estructura real de columnsResult
            // console.log(columnsResult);

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
        console.error('❌ Error en getSQLPrompt:', err);
        throw err;
    } finally {
        if (connection) await connection.end();
    }
}

// Ejecución principal para Flowise
try {
    await getSQLPrompt();
    return sqlSchemaPrompt;
} catch (error) {
    return `❌ Error generando listado de tablas: ${error.message}`;
}
