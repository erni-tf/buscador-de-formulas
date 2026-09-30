"""Genera el archivo único `dist/Buscador_de_Formulas.html`.

Lee la base de datos, dibuja todas las fórmulas (esquelética y
semidesarrollada), y las empaqueta junto con el buscador, el parser y la
librería de dibujo en un solo HTML que funciona sin internet.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(RAIZ / "src"))

import catalogo  # noqa: E402
from dibujo_condensada import svg_condensada  # noqa: E402
from dibujo_esqueleto import svg_condensado_anillo, svg_esqueleto  # noqa: E402
from rdkit import Chem, RDLogger  # noqa: E402
from rdkit.Chem import rdMolDescriptors  # noqa: E402

RDLogger.DisableLog("rdApp.*")

PLANTILLA = RAIZ / "src" / "plantilla"
DIST = RAIZ / "dist"
SALIDA = DIST / "Buscador_de_Formulas.html"

# Las fichas genéricas (cetona R-CO-R'...) llevan etiquetas R en vez de línea
# ondulada; los sustituyentes (vinilo, fenilo...) llevan el punto de unión.
FAMILIAS_CON_R = {"grupos funcionales", "conceptos"}


def formula_molecular(smiles: str) -> str:
    """Fórmula molecular. Vacía para las fichas genéricas (llevan comodines)."""
    if "*" in smiles:
        return ""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""
    return rdMolDescriptors.CalcMolFormula(mol)


def es_anillo(smiles: str) -> bool:
    mol = Chem.MolFromSmiles(smiles)
    return bool(mol and mol.GetRingInfo().NumRings() > 0)


def tiene_nitrogeno(smiles: str) -> bool:
    """True si hay algún nitrógeno con par libre (los dos puntitos)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    return any(
        a.GetAtomicNum() == 7
        and a.GetFormalCharge() == 0
        and all(b.GetBondType() == Chem.BondType.SINGLE for b in a.GetBonds())
        for a in mol.GetAtoms()
    )


def tiene_aromatico(smiles: str) -> bool:
    """True si hay algún anillo aromático (se puede dibujar con círculo)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    return any(
        all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in anillo)
        for anillo in mol.GetRingInfo().AtomRings()
    )


def dibujar(entrada: catalogo.Entrada) -> dict[str, str]:
    """Devuelve {'esq': svg, 'con': svg, ...} para una entrada.

    Si las dos vistas salen idénticas (benceno, un sustituyente sencillo…),
    se quita la segunda para no repetir el mismo dibujo. Si hay anillos
    aromáticos se generan además las versiones con círculo ('-circ'). La
    fórmula semidesarrollada lleva variantes para los conmutadores: '-cooh'
    (los COOH separados en C(=O)-OH), '-ch2' (los CH2 seguidos juntos) y
    '-ch3' (los CH3 repetidos del mismo átomo en un grupo).
    """
    imagenes: dict[str, str] = {}

    con_r = any(f in FAMILIAS_CON_R for f in entrada.familias)
    imagenes["esq"] = svg_esqueleto(entrada.smiles, etiquetas_dummy=con_r)

    aromatico = tiene_aromatico(entrada.smiles)
    if aromatico:
        imagenes["esq-circ"] = svg_esqueleto(
            entrada.smiles, etiquetas_dummy=con_r, circulos=True
        )

    if entrada.condensada:
        for clave, cooh, ch2, ch3 in (
            ("con", False, False, False),
            ("con-cooh", True, False, False),
            ("con-ch2", False, True, False),
            ("con-cooh-ch2", True, True, False),
            ("con-ch3", False, False, True),
            ("con-cooh-ch3", True, False, True),
            ("con-ch2-ch3", False, True, True),
            ("con-cooh-ch2-ch3", True, True, True),
        ):
            imagenes[clave] = svg_condensada(
                entrada.condensada, cooh_separado=cooh, ch2_agrupados=ch2,
                ch3_agrupados=ch3,
            )
    elif es_anillo(entrada.smiles):
        for clave, cooh, ch3 in (
            ("con", False, False),
            ("con-cooh", True, False),
            ("con-ch3", False, True),
            ("con-cooh-ch3", True, True),
        ):
            imagenes[clave] = svg_condensado_anillo(
                entrada.smiles, cooh_separado=cooh, ch3_agrupados=ch3
            )
            if aromatico:
                imagenes[f"con-circ{clave[len('con'):]}"] = svg_condensado_anillo(
                    entrada.smiles, circulos=True, cooh_separado=cooh,
                    ch3_agrupados=ch3,
                )

    return _quitar_repetidas(imagenes)


def _quitar_repetidas(imagenes: dict[str, str]) -> dict[str, str]:
    """Deja una sola copia de cada dibujo, con el nombre más simple.

    Se recorre en orden (esq, con, variantes, circulares): si un dibujo ya
    había salido, se descarta. Así el benceno no repite la vista «con» y una
    fórmula sin COOH no lleva la variante «con-cooh» (sería idéntica).
    """
    vistos: set[str] = set()
    salida: dict[str, str] = {}
    for clave, svg in imagenes.items():
        if not svg or svg in vistos:
            continue
        vistos.add(svg)
        salida[clave] = svg
    return salida


# ---------------------------------------------------------------------------
#  Las vistas: qué dibujo enseña cada combinación de los conmutadores
#
#  Los tres chips que cambian de dibujo son el anillo, los COOH y los CH2
#  seguidos (el de los 2 puntos de la N no elige dibujo: sólo enseña o
#  esconde un grupo dentro del SVG). De aquí sale una tabla por entrada, con
#  una casilla por combinación, para que la aplicación no tenga que saber
#  nada del espacio de claves: lee la que le toca.
# ---------------------------------------------------------------------------

ANILLOS = ("rayas", "circulo")
COOHS = ("junto", "separado")
CH2S = ("sueltos", "agrupados")
CH3S = ("sueltos", "agrupados")

# Lo que se prueba en cada estado de los chips de la fórmula: primero la
# variante más completa y después las que van soltando conmutadores; el COOH
# manda sobre el CH2 y el CH2 sobre el CH3, como cuando lo decidía la
# aplicación. Las claves van en el orden del nombre de la variante
# («con-cooh-ch2-ch3»).
_SUFIJOS = {
    ("junto", "sueltos", "sueltos"): ("",),
    ("junto", "sueltos", "agrupados"): ("ch3", ""),
    ("junto", "agrupados", "sueltos"): ("ch2", ""),
    ("junto", "agrupados", "agrupados"): ("ch2-ch3", "ch2", "ch3", ""),
    ("separado", "sueltos", "sueltos"): ("cooh", ""),
    ("separado", "sueltos", "agrupados"): ("cooh-ch3", "cooh", "ch3", ""),
    ("separado", "agrupados", "sueltos"): ("cooh-ch2", "cooh", "ch2", ""),
    ("separado", "agrupados", "agrupados"): (
        "cooh-ch2-ch3",
        "cooh-ch2",
        "cooh-ch3",
        "cooh",
        "ch2-ch3",
        "ch2",
        "ch3",
        "",
    ),
}


def resolver_vistas(imagenes: dict[str, str]) -> dict[str, dict[str, str]]:
    """Qué dibujo enseña cada combinación de los conmutadores.

    Devuelve `estructura` (rayas o círculo) y `semidesarrollada` (las
    dieciséis combinaciones de anillo, COOH, CH2 y CH3), que no aparece si la
    entrada no tiene ninguna variante `con*`. La política es la de siempre:
    primero el anillo —el círculo, si el conmutador está encendido y existe—
    y después la variante de la fórmula más completa que se haya dibujado, de
    modo que un conmutador sin variante no cambia el dibujo y toda casilla
    tiene una clave que existe.
    """
    vistas: dict[str, dict[str, str]] = {
        "estructura": {
            "rayas": "esq",
            "circulo": "esq-circ" if "esq-circ" in imagenes else "esq",
        }
    }
    if "con" not in imagenes:
        return vistas

    semidesarrollada: dict[str, str] = {}
    for anillo in ANILLOS:
        circulares = ("-circ", "") if anillo == "circulo" else ("",)
        for cooh in COOHS:
            for ch2 in CH2S:
                for ch3 in CH3S:
                    clave = ""
                    for circular in circulares:
                        for sufijo in _SUFIJOS[(cooh, ch2, ch3)]:
                            candidata = "con" + circular + (f"-{sufijo}" if sufijo else "")
                            if candidata in imagenes:
                                clave = candidata
                                break
                        if clave:
                            break
                    semidesarrollada[f"{anillo}-{cooh}-{ch2}-{ch3}"] = clave
    vistas["semidesarrollada"] = semidesarrollada
    return vistas


def conmutadores(vistas: dict[str, dict[str, str]]) -> dict[str, bool]:
    """Qué conmutadores de la fórmula cambian el dibujo (para sus chips)."""
    claves = list(vistas.get("semidesarrollada", {}).values())
    return {
        "cooh": any("-cooh" in clave for clave in claves),
        "ch2": any("-ch2" in clave for clave in claves),
        "ch3": any("-ch3" in clave for clave in claves),
    }


def main() -> int:
    entradas = catalogo.cargar()
    problemas = catalogo.validar(entradas)
    if problemas:
        print("La base de datos tiene problemas:")
        for p in problemas:
            print("  -", p)
        return 1

    bloques_imagenes: list[str] = []
    datos: list[dict] = []

    for e in entradas:
        imagenes = dibujar(e)
        vistas = resolver_vistas(imagenes)
        for tipo, svg in imagenes.items():
            bloques_imagenes.append(
                f'<div class="grafico" id="img-{tipo}-{e.id}" hidden>{svg}</div>'
            )
        datos.append({
            "id": e.id,
            "nombre": e.nombre,
            "nombres": e.todos_los_nombres,
            "alias": e.alias,
            "familias": e.familias,
            "tipo": e.tipo,
            "formula": formula_molecular(e.smiles),
            "notas": e.notas,
            "condensada_texto": e.condensada,
            "vistas": vistas,
            "circulos": tiene_aromatico(e.smiles),
            "puntos": tiene_nitrogeno(e.smiles),
            **conmutadores(vistas),
        })

    familias = [f for f in catalogo.FAMILIAS if any(f in d["familias"] for d in datos)]

    html = (PLANTILLA / "app.html").read_text(encoding="utf-8")
    html = html.replace("/*__CSS__*/", (PLANTILLA / "app.css").read_text(encoding="utf-8"))
    html = html.replace("__IMAGENES__", "\n".join(bloques_imagenes))
    html = html.replace(
        "/*__SMILESDRAWER__*/",
        (PLANTILLA / "smiles-drawer.min.js").read_text(encoding="utf-8"),
    )
    html = html.replace("/*__PARSER__*/", (PLANTILLA / "parser.js").read_text(encoding="utf-8"))
    html = html.replace("/*__APP__*/", (PLANTILLA / "app.js").read_text(encoding="utf-8"))
    html = html.replace(
        "__DATOS__",
        json.dumps(
            {"entradas": datos, "familias": familias},
            ensure_ascii=False,
            separators=(",", ":"),
        ),
    )

    DIST.mkdir(exist_ok=True)
    SALIDA.write_text(html, encoding="utf-8")

    kb = SALIDA.stat().st_size / 1024
    print(f"Generado: {SALIDA}")
    print(f"  entradas: {len(datos)}   imágenes: {len(bloques_imagenes)}")
    print(f"  tamaño: {kb:,.0f} KB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
