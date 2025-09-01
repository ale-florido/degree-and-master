FROM node:20
# Usa la imagen oficial de Node.js versión 20 como base.

WORKDIR /usr/src/app
# Establece el directorio de trabajo dentro del contenedor.

COPY package*.json ./
# Copia los archivos package.json y package-lock.json (si existe) al directorio de trabajo.

RUN npm install
# Instala las dependencias definidas en package.json.

COPY . .
# Copia el resto de los archivos de tu proyecto al directorio de trabajo en el contenedor.

EXPOSE 3000
# Expone el puerto 3000 para que el contenedor pueda recibir conexiones externas en ese puerto.

CMD ["node", "server.js"]
# Comando por defecto al iniciar el contenedor: ejecuta server.js con Node.js.
