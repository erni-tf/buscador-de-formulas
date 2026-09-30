"""Comprueba la base de datos, que todas las entradas se dibujan bien y que
la tabla de vistas (qué dibujo enseña cada combinación de conmutadores) es
coherente con lo dibujado.

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


# ---------------------------------------------------------------------------
#  La tabla de vistas
#
#  Cada caso le da al resolvedor un diccionario de claves dibujadas (como el
#  que sale de `construir.dibujar`) y dice qué casilla tiene que elegir. Son
#  los casos que fijan la política: el anillo manda sobre la fórmula, la
#  variante más completa manda, y un conmutador sin variante no cambia nada.
# ---------------------------------------------------------------------------

# (nombre, claves dibujadas, casilla del anillo, clave esperada)
CASOS_ESTRUCTURA = [
    ("con círculo dibujado", ["esq", "esq-circ"], "circulo", "esq-circ"),
    ("sin círculo dibujado", ["esq"], "circulo", "esq"),
    ("las rayas siempre están", ["esq", "esq-circ"], "rayas", "esq"),
]

# (nombre, claves dibujadas, combinación, clave esperada)
CASOS_FORMULA = [
    ("sin variantes dibujadas", ["esq", "con"], "rayas-junto-sueltos", "con"),
    ("el COOH pedido manda", ["esq", "con", "con-cooh"], "rayas-separado-sueltos", "con-cooh"),
    ("el COOH no dibujado cae al con", ["esq", "con"], "rayas-separado-sueltos", "con"),
    ("el CH2 pedido manda", ["esq", "con", "con-ch2"], "rayas-junto-agrupados", "con-ch2"),
    ("el CH2 no dibujado cae al con", ["esq", "con", "con-cooh"], "rayas-junto-agrupados", "con"),
    (
        "la variante completa manda",
        ["esq", "con", "con-cooh", "con-ch2", "con-cooh-ch2"],
        "rayas-separado-agrupados",
        "con-cooh-ch2",
    ),
    (
        "el COOH manda sobre el CH2",
        ["esq", "con", "con-cooh", "con-ch2"],
        "rayas-separado-agrupados",
        "con-cooh",
    ),
    (
        "el círculo manda sobre la fórmula",
        ["esq", "esq-circ", "con", "con-circ", "con-cooh"],
        "circulo-separado-sueltos",
        "con-circ",
    ),
    (
        "el círculo con COOH, si existe",
        ["esq", "esq-circ", "con", "con-circ", "con-circ-cooh", "con-cooh"],
        "circulo-separado-sueltos",
        "con-circ-cooh",
    ),
    (
        "el círculo sin variante de CH2",
        ["esq", "esq-circ", "con", "con-circ", "con-circ-cooh"],
        "circulo-junto-agrupados",
        "con-circ",
    ),
    (
        "sin círculo dibujado, la fórmula sigue",
        ["esq", "con", "con-cooh"],
        "circulo-separado-sueltos",
        "con-cooh",
    ),
]


def _revisar_vistas(entradas: list, dibujos: dict[str, dict[str, str]]) -> list[str]:
    """La política del resolvedor y la tabla de todas las entradas."""
    fallos: list[str] = []

    for nombre, claves, anillo, esperada in CASOS_ESTRUCTURA:
        vistas = construir.resolver_vistas({clave: "x" for clave in claves})
        obtenida = vistas["estructura"][anillo]
        if obtenida != esperada:
            fallos.append(f"{nombre}: se esperaba {esperada!r}, sale {obtenida!r}")

    for nombre, claves, combinacion, esperada in CASOS_FORMULA:
        vistas = construir.resolver_vistas({clave: "x" for clave in claves})
        obtenida = vistas["semidesarrollada"][combinacion]
        if obtenida != esperada:
            fallos.append(f"{nombre}: se esperaba {esperada!r}, sale {obtenida!r}")

    # sin ninguna variante «con» no hay sub-tabla de fórmula
    if "semidesarrollada" in construir.resolver_vistas({"esq": "x"}):
        fallos.append("sin variantes «con»: no debería haber sub-tabla")

    # y el catálogo entero: cada casilla apunta a un dibujo que existe
    for e in entradas:
        imagenes = dibujos.get(e.id)
        if imagenes is None:
            continue
        vistas = construir.resolver_vistas(imagenes)
        for clave in vistas["estructura"].values():
            if clave not in imagenes:
                fallos.append(f"{e.id}: la estructura apunta a {clave!r}, que no existe")
        if ("semidesarrollada" in vistas) != ("con" in imagenes):
            fallos.append(f"{e.id}: la sub-tabla no concuerda con las variantes dibujadas")
        for clave in vistas.get("semidesarrollada", {}).values():
            if clave not in imagenes:
                fallos.append(f"{e.id}: la fórmula apunta a {clave!r}, que no existe")
        esperados = {
            "cooh": "con-cooh" in imagenes or "con-circ-cooh" in imagenes,
            "ch2": "con-ch2" in imagenes or "con-cooh-ch2" in imagenes,
        }
        if construir.conmutadores(vistas) != esperados:
            fallos.append(f"{e.id}: los chips no concuerdan con los dibujos")

    return fallos


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
    dibujos: dict[str, dict[str, str]] = {}
    for e in entradas:
        try:
            imagenes = construir.dibujar(e)
        except Exception as exc:  # noqa: BLE001
            fallos.append(f"{e.id}: {exc}")
            continue
        dibujos[e.id] = imagenes
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

    fallos_vistas = _revisar_vistas(entradas, dibujos)

    if fallos or fallos_vistas:
        if fallos:
            print(f"\nFALLOS DE DIBUJO ({len(fallos)}):")
            for f in fallos:
                print("  -", f)
        if fallos_vistas:
            print(f"\nFALLOS DE VISTAS ({len(fallos_vistas)}):")
            for f in fallos_vistas:
                print("  -", f)
        return 1

    # Comprobación extra: todas las fórmulas moleculares se pueden calcular
    for e in entradas:
        if "*" in e.smiles:
            continue
        if Chem.MolFromSmiles(e.smiles) is None:
            print(f"\nSMILES inválido en {e.id}")
            return 1

    print(f"\nVistas: {len(CASOS_ESTRUCTURA) + len(CASOS_FORMULA)} casos de política "
          f"y {len(dibujos)} entradas revisadas.")
    print("Base de datos, dibujos y vistas correctos.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
