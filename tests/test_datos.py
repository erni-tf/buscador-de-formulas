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
from dibujo_condensada import _agrupar_ch3, analizar  # noqa: E402
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
    ("sin variantes dibujadas", ["esq", "con"], "rayas-junto-sueltos-sueltos", "con"),
    ("el COOH pedido manda", ["esq", "con", "con-cooh"], "rayas-separado-sueltos-sueltos", "con-cooh"),
    ("el COOH no dibujado cae al con", ["esq", "con"], "rayas-separado-sueltos-sueltos", "con"),
    ("el CH2 pedido manda", ["esq", "con", "con-ch2"], "rayas-junto-agrupados-sueltos", "con-ch2"),
    ("el CH2 no dibujado cae al con", ["esq", "con", "con-cooh"], "rayas-junto-agrupados-sueltos", "con"),
    ("el CH3 pedido manda", ["esq", "con", "con-ch3"], "rayas-junto-sueltos-agrupados", "con-ch3"),
    ("el CH3 no dibujado cae al con", ["esq", "con", "con-cooh"], "rayas-junto-sueltos-agrupados", "con"),
    (
        "la variante completa manda",
        ["esq", "con", "con-cooh", "con-ch2", "con-cooh-ch2", "con-ch3",
         "con-cooh-ch3", "con-ch2-ch3", "con-cooh-ch2-ch3"],
        "rayas-separado-agrupados-agrupados",
        "con-cooh-ch2-ch3",
    ),
    (
        "el COOH manda sobre el CH2 y el CH3",
        ["esq", "con", "con-cooh", "con-ch2", "con-ch3"],
        "rayas-separado-agrupados-agrupados",
        "con-cooh",
    ),
    (
        "el CH2 manda sobre el CH3",
        ["esq", "con", "con-ch2", "con-ch3"],
        "rayas-junto-agrupados-agrupados",
        "con-ch2",
    ),
    (
        "el círculo manda sobre la fórmula",
        ["esq", "esq-circ", "con", "con-circ", "con-cooh"],
        "circulo-separado-sueltos-sueltos",
        "con-circ",
    ),
    (
        "el círculo con COOH, si existe",
        ["esq", "esq-circ", "con", "con-circ", "con-circ-cooh", "con-cooh"],
        "circulo-separado-sueltos-sueltos",
        "con-circ-cooh",
    ),
    (
        "el círculo con CH3, si existe",
        ["esq", "esq-circ", "con", "con-circ", "con-ch3", "con-circ-ch3"],
        "circulo-junto-sueltos-agrupados",
        "con-circ-ch3",
    ),
    (
        "el círculo sin variante de CH2",
        ["esq", "esq-circ", "con", "con-circ", "con-circ-cooh"],
        "circulo-junto-agrupados-sueltos",
        "con-circ",
    ),
    (
        "sin círculo dibujado, la fórmula sigue",
        ["esq", "con", "con-cooh"],
        "circulo-separado-sueltos-sueltos",
        "con-cooh",
    ),
]

# ---------------------------------------------------------------------------
#  El conmutador «CH3»: el agrupado que hace dibujo_condensada._agrupar_ch3
#
#  Cada caso da la fórmula del DSL y el resultado esperado (el mismo DSL,
#  como lo devolvería el serializador de abajo).
# ---------------------------------------------------------------------------

CASOS_CH3 = [
    ("isopropilo", "CH3-CH[CH3]-O-C[=O]-H", "(CH3)2CH-O-C[=O]-H"),
    ("dimetilformamida", "H-C[=O]-N[CH3]-CH3", "H-C[=O]-N(CH3)2"),
    ("tert-butilamina", "CH3-C{NH2}[CH3]-CH3", "(CH3)3C-NH2"),
    ("dietilmetilamina", "CH3-CH2-N[CH3]-CH2-CH3", "(CH3CH2)2N-CH3"),
    ("éter dietílico", "CH3-CH2-O-CH2-CH3", "(CH3CH2)2O"),
    ("ácido ramificado", "CH3-CH[CH3]-CH2-COOH", "(CH3)2CH-CH2-COOH"),
    ("alqueno", "CH2=CH-CH[CH3]-CH3", "CH2=CH-CH(CH3)2"),
    ("tres metilos", "CH3-C{CH3}[CH3]-CH{OH}-CH2-Br", "(CH3)3C-CH{OH}-CH2-Br"),
    ("una cadena seguida no se agrupa", "CH3-CH2-CH2-CH3", "CH3-CH2-CH2-CH3"),
]


def _texto_cadena(cadena) -> str:
    """Vuelve a escribir la cadena en el DSL (para los casos del CH3)."""
    partes = []
    for nodo, enlace in cadena:
        texto = nodo.grupo
        for lado, etiqueta in ((nodo.arriba, ("{", "}")), (nodo.abajo, ("[", "]"))):
            for rama in lado:
                contenido = _texto_cadena(rama.cadena)
                if rama.enlace != "-":
                    contenido = rama.enlace + contenido
                texto += etiqueta[0] + contenido + etiqueta[1]
        partes.append(texto + enlace)
    return "".join(partes)


# Claves que puede generar construir.dibujar() para la fórmula
# semidesarrollada: cadena (con ch3) o anillo (sin ch2).
CLAVES_CONOCIDAS = {
    "con", "con-cooh", "con-ch2", "con-cooh-ch2", "con-ch3", "con-cooh-ch3",
    "con-ch2-ch3", "con-cooh-ch2-ch3",
    "con-circ", "con-circ-cooh", "con-circ-ch3", "con-circ-cooh-ch3",
}


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

    # el agrupado de los CH3 (y de los etilos) de las fórmulas de cadena
    for nombre, formula, esperada in CASOS_CH3:
        obtenida = _texto_cadena(_agrupar_ch3(analizar(formula.strip("-").strip())))
        if obtenida != esperada:
            fallos.append(f"CH3 ({nombre}): se esperaba {esperada!r}, sale {obtenida!r}")

    # sin ninguna variante «con» no hay sub-tabla de fórmula
    if "semidesarrollada" in construir.resolver_vistas({"esq": "x"}):
        fallos.append("sin variantes «con»: no debería haber sub-tabla")

    # y el catálogo entero: cada casilla apunta a un dibujo que existe
    for e in entradas:
        imagenes = dibujos.get(e.id)
        if imagenes is None:
            continue
        for clave in imagenes:
            if clave not in ("esq", "esq-circ") and clave not in CLAVES_CONOCIDAS:
                fallos.append(f"{e.id}: clave dibujada desconocida {clave!r}")
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
            "ch3": "con-ch3" in imagenes or "con-circ-ch3" in imagenes,
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

    print(f"\nVistas: {len(CASOS_ESTRUCTURA) + len(CASOS_FORMULA)} casos de política, "
          f"{len(CASOS_CH3)} de agrupado del CH3 y {len(dibujos)} entradas revisadas.")
    print("Base de datos, dibujos y vistas correctos.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
