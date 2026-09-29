"""Comprueba la base de datos y que todas las entradas se dibujan bien.

Se ejecuta con:  .venv/bin/python tests/test_datos.py
"""

from __future__ import annotations

import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(RAIZ / "src"))

import catalogo  # noqa: E402
import construir  # noqa: E402
from rdkit import Chem, RDLogger  # noqa: E402

RDLogger.DisableLog("rdApp.*")


def main() -> int:
    entradas = catalogo.cargar()
    problemas = catalogo.validar(entradas)

    print(f"Entradas: {len(entradas)}")
    print(f"  compuestos: {sum(1 for e in entradas if e.tipo == 'compuesto')}")
    print(f"  fichas:     {sum(1 for e in entradas if e.tipo == 'ficha')}")

    if problemas:
        print(f"\nPROBLEMAS EN LA BASE ({len(problemas)}):")
        for p in problemas:
            print("  -", p)
        return 1

    # Dibujar todas las entradas (detecta errores de renderizado)
    fallos: list[str] = []
    sin_condensada: list[str] = []
    for e in entradas:
        try:
            imagenes = construir.dibujar(e)
        except Exception as exc:  # noqa: BLE001
            fallos.append(f"{e.id}: {exc}")
            continue
        if not imagenes.get("esq"):
            fallos.append(f"{e.id}: sin vista esquelética")
        if not imagenes.get("con"):
            sin_condensada.append(e.id)
        for tipo, svg in imagenes.items():
            if svg and "<svg" not in svg:
                fallos.append(f"{e.id} ({tipo}): el SVG no parece válido")

    if sin_condensada:
        print(f"\nSin vista semidesarrollada ({len(sin_condensada)}) — "
              f"se usará la esquelética:")
        for i in sin_condensada:
            print("  ·", i)

    if fallos:
        print(f"\nFALLOS DE DIBUJO ({len(fallos)}):")
        for f in fallos:
            print("  -", f)
        return 1

    # Comprobación extra: todas las fórmulas moleculares se pueden calcular
    for e in entradas:
        if "*" in e.smiles:
            continue
        if Chem.MolFromSmiles(e.smiles) is None:
            print(f"\nSMILES inválido en {e.id}")
            return 1

    print("\nBase de datos y dibujos correctos.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
