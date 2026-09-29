"""Comprueba que el parser (JavaScript) reproduce los SMILES del catálogo.

Ejecuta `parser.js` con Node y compara cada resultado con el SMILES curado,
usando RDKit como árbitro (canonicaliza las dos cadenas antes de comparar).
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(RAIZ / "src"))

import catalogo  # noqa: E402
from rdkit import Chem, RDLogger  # noqa: E402

RDLogger.DisableLog("rdApp.*")

PARSER = RAIZ / "src" / "plantilla" / "parser.js"

# Nombres de los apuntes a los que les falta un localizador. El parser hace una
# suposición (y avisa de ella), así que pueden dar un isómero distinto. El
# catálogo sí tiene la respuesta correcta para estos.
AMBIGUOS = {
    "6-bromo-2-metilheptino",  # falta el localizador del triple enlace
    "feniletanol",             # falta el localizador del fenilo
    "1-bromo-3,3-dimetilbutanol",  # falta el localizador del OH
}

GUION = """
const P = require(process.env.PARSER_JS);
const casos = JSON.parse(require('fs').readFileSync(0, 'utf8'));
const salida = casos.map(function (c) {
  var r = P.construir(c.nombre);
  return { id: c.id, nombre: c.nombre, esperado: c.esperado, resultado: r };
});
process.stdout.write(JSON.stringify(salida));
"""


def canon(smiles: str) -> str | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol)


def node_bin() -> str:
    return shutil.which("node") or "/home/ernesttf/.local/share/pi-node/node-v22.23.2-linux-x64/bin/node"


def main() -> int:
    entradas = [e for e in catalogo.cargar() if e.tipo == "compuesto"]
    casos = [{"id": e.id, "nombre": e.nombre, "esperado": e.smiles} for e in entradas]
    # Además, todas las variantes (alias) que sean nombres sistemáticos
    for e in entradas:
        for alias in e.alias:
            casos.append({"id": e.id, "nombre": alias, "esperado": e.smiles})

    proc = subprocess.run(
        [node_bin(), "-e", GUION],
        input=json.dumps(casos),
        capture_output=True,
        text=True,
        env={**os.environ, "PARSER_JS": str(PARSER)},
    )
    if proc.returncode != 0:
        print("Error ejecutando Node:")
        print(proc.stderr)
        return 2

    resultados = json.loads(proc.stdout)
    fallos: list[tuple[str, str, str, str]] = []
    supuestos: list[tuple[str, str, str]] = []
    for r in resultados:
        esperado = canon(r["esperado"])
        res = r["resultado"]
        if not res.get("ok"):
            fallos.append((r["nombre"], r["esperado"], "—", res.get("error", "?")))
            continue
        obtenido = canon(res["smiles"])
        if obtenido != esperado:
            if r["nombre"] in AMBIGUOS:
                supuestos.append((r["nombre"], esperado or "", obtenido or res["smiles"]))
            else:
                fallos.append((r["nombre"], esperado, obtenido or res["smiles"], "no coincide"))

    total = len(resultados)
    print(f"Casos: {total}   correctos: {total - len(fallos) - len(supuestos)}   "
          f"fallos: {len(fallos)}   suposiciones documentadas: {len(supuestos)}")
    if supuestos:
        print()
        print("Nombres incompletos de los apuntes (el parser avisa de la suposición):")
        for nombre, esp, obt in supuestos:
            print(f"  · {nombre}:  isómero asumido {obt}  (el correcto es {esp})")
    if fallos:
        print()
        for nombre, esp, obt, motivo in fallos:
            print(f"  ✗ {nombre}")
            print(f"      esperado: {esp}")
            print(f"      obtenido: {obt}   ({motivo})")
        return 1
    print("Todos los nombres sistemáticos del catálogo se interpretan correctamente.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
