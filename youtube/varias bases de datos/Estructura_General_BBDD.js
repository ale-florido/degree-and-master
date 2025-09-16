const mysql = require('mysql2/promise');

// Configuración base común
const COMMON_CONFIG = {
    host: 'mysql-sakila-test', //CAMBIA POR EL NOMBRE DEL CONTENEDOR CON LAS BASES DE DATOS
    user: 'root',
    password: 'password',
    port: 3306,
    timeout: 5000
};

// Lista de bases de datos disponibles
const AVAILABLE_DATABASES = ['world', 'sakila'];

async function getDatabaseSchemaPrompt(databaseName) {
    // Validar que la base de datos esté en la lista disponible
    if (!AVAILABLE_DATABASES.includes(databaseName)) {
        return `❌ La base de datos '${databaseName}' no está disponible. Bases de datos disponibles: ${AVAILABLE_DATABASES.join(', ')}`;
    }

    let connection;
    
    try {
        connection = await mysql.createConnection({
            ...COMMON_CONFIG,
            database: databaseName
        });

        // Obtener todos los nombres de tablas
        const [tablesResult] = await connection.query(`
            SELECT table_name
            FROM information_schema.tables
            WHERE table_schema = ? AND table_type = 'BASE TABLE'
            ORDER BY table_name ASC
        `, [databaseName]);

        if (tablesResult.length === 0) {
            return `⚠️ No se encontraron tablas en la base de datos \`${databaseName}\`.`;
        }

        // Construir prompt con nombres de tablas en Markdown
        let prompt = `✅ **Lista de tablas en la base de datos \`${databaseName}\`:**\n\n`;

        for (const tableRow of tablesResult) {
            const tableName = tableRow.table_name || tableRow.TABLE_NAME;
            prompt += `- \`${tableName}\`\n`;
        }

        return prompt;

    } catch (err) {
        console.error(`❌ Error conectando a la base de datos '${databaseName}':`, err.message);
        throw new Error(`Error al conectar con la base de datos '${databaseName}': ${err.message}`);
    } finally {
        if (connection) await connection.end();
    }
}

async function getMultipleDatabaseSchemas(databaseNames) {
    const results = {};
    
    for (const dbName of databaseNames) {
        try {
            results[dbName] = await getDatabaseSchemaPrompt(dbName);
        } catch (error) {
            results[dbName] = `❌ Error: ${error.message}`;
        }
    }
    
    return results;
}

// Ejecución para Flowise
try {
    const multipleResults = await getMultipleDatabaseSchemas(['world', 'sakila']);
    return JSON.stringify(multipleResults, null, 2);

} catch (error) {
    return `❌ Error generando listado de tablas: ${error.message}`;
}
