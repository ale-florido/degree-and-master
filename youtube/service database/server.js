const express = require('express');
// Importa el framework Express.

const { exec } = require('child_process');
// Importa la función exec para ejecutar comandos del sistema.

const app = express();
// Crea una instancia de la aplicación Express.

app.get('/run-mysql-script', (req, res) => {
  // Define una ruta GET en /run-mysql-script.
  exec('node mysql-script.js', (error, stdout, stderr) => {
    // Ejecuta el script mysql-script.js como un proceso hijo.
    if (error) {
      // Si ocurre un error al ejecutar el script, responde con error 500.
      res.status(500).json({ error: error.message });
      return;
    }
    try {
      // Si no hay error, envía la salida estándar (stdout) como respuesta.
      res.send(stdout);
    } catch (e) {
      // Si ocurre un error al procesar la salida, responde con error 500.
      res.status(500).json({ error: 'Error procesando la salida del script', details: stdout });
    }
  });
});

app.listen(3000, () => {
  // Inicia el servidor Express en el puerto 3000.
  console.log('Servidor Express escuchando en puerto 3000');
});

