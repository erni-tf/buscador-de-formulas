#!/usr/bin/env bash
# ---------------------------------------------------------------------------
#  Regenera el archivo único para Laura a partir de datos/compuestos.toml
#  Uso:  ./generar.sh
# ---------------------------------------------------------------------------
set -e
cd "$(dirname "$0")"

# Usa el entorno virtual del proyecto si existe; si no, el python del sistema.
PY=".venv/bin/python"
if [ ! -x "$PY" ]; then
  PY="python3"
fi

echo "1/4  Comprobando la base de datos…"
"$PY" tests/test_datos.py

echo
echo "2/4  Comprobando la postura de las moléculas…"
"$PY" tests/test_postura.py

echo
echo "3/4  Comprobando el parser de nombres…"
"$PY" tests/test_parser.py

echo
echo "4/4  Generando el archivo…"
"$PY" src/construir.py

echo
echo "Listo. El archivo que hay que enviar a Laura es:"
echo "  dist/Buscador_de_Formulas.html"
