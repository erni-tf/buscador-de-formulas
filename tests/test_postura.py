"""Comprueba la postura: qué regla reclama cada molécula y cómo la deja.

Dos preguntas, dos listas:

  * la **tabla de reglas**: cada compuesto de los apuntes tiene que ser
    reclamado por su regla (y no por otra, ni por ninguna);
  * la **geometría**: las posturas con invariante claro se comprueban
    midiendo el conformador (el O del pentenol debajo de la cadena, el del
    anisol a la derecha del anillo, los dos anillos del éster a los lados).

Se ejecuta con:  .venv/bin/python tests/test_postura.py
"""

from __future__ import annotations

import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(RAIZ / "src"))

import dibujo_esqueleto  # noqa: E402
import postura  # noqa: E402
from rdkit import Chem, RDLogger  # noqa: E402
from rdkit.Chem.Draw import rdMolDraw2D  # noqa: E402

RDLogger.DisableLog("rdApp.*")


def _colocada(smiles: str, con_grupos: bool = False):
    """La molécula preparada como la deja el dibujo, ya colocada.

    Se prepara con CoordGen, igual que `_dibujar`. Con `con_grupos` pasa
    antes por la tabla de grupos de la vista condensada —el mismo ayudante
    `_colapsar_grupos` que usa el dibujo—, que es lo que hace que el éster
    diarílico o las anilinas lleven rótulo.
    """
    mol = rdMolDraw2D.PrepareMolForDrawing(Chem.MolFromSmiles(smiles))
    etiquetas: dict[int, str] = {}
    anclas: dict[int, int | str] = {}
    if con_grupos:
        mol, etiquetas, anclas = dibujo_esqueleto._colapsar_grupos(mol)
        mol = rdMolDraw2D.PrepareMolForDrawing(mol)
    registro = postura.postura(mol, etiquetas, anclas)
    return mol, etiquetas, anclas, registro


def _unico(mol: Chem.Mol, numero: int) -> int:
    """Índice del único átomo de ese elemento en la molécula."""
    return [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == numero][0]


# (nombre, SMILES, ¿con la tabla de grupos?, regla esperada)
CASOS = [
    # las reglas de compuesto
    ("2-penten-3-ol", "CCC(O)=CC", False, "penten_3_ol"),
    ("metoxibenceno", "COc1ccccc1", False, "metoxibenceno"),
    ("1-etil-2-metilciclobutano", "CCC1CCC1C", False, "ciclobutano"),
    ("1-etil-2-metilciclobutano (etilo escrito)", "CCC1CCC1C", True, "ciclobutano"),
    ("ciclohexanocarbaldehido", "O=CC1CCCCC1", False, "ciclohexanocarbaldehido"),
    ("acido-2-ciclohexilpropanoico", "CC(C(=O)O)C1CCCCC1", False, "acido_ciclohexilpropanoico"),
    ("acido-2-ciclohexilpropanoico (rótulo)", "CC(C(=O)O)C1CCCCC1", True, "acido_ciclohexilpropanoico"),
    ("benzoato-de-fenilo", "O=C(Oc1ccccc1)c1ccccc1", True, "ester_diarilico"),
    # el grupo «sustituyentes»: los apuntes los dibujan hacia la derecha
    ("etinilbenceno", "C#Cc1ccccc1", False, "sustituyentes"),
    ("1-fenil-2-butino", "CC#CCc1ccccc1", False, "sustituyentes"),
    ("2-feniletanol", "OCCc1ccccc1", False, "sustituyentes"),
    ("bencilamina", "NCc1ccccc1", False, "sustituyentes"),
    ("bencenoato", "[O-]C(=O)c1ccccc1", True, "sustituyentes"),
    ("N-etilanilina", "CCNc1ccccc1", True, "sustituyentes"),
    ("N,N-dietilanilina", "CCN(CC)c1ccccc1", True, "sustituyentes"),
]


def _revisar_pentenol() -> str | None:
    """El -OH del 2-penten-3-ol queda debajo de la cadena."""
    mol, _, _, _ = _colocada("CCC(O)=CC")
    conf = mol.GetConformer()
    o = _unico(mol, 8)
    c3 = mol.GetAtomWithIdx(o).GetNeighbors()[0].GetIdx()
    if conf.GetAtomPosition(o).y >= conf.GetAtomPosition(c3).y - 1.0:
        return "el O del 2-penten-3-ol no queda debajo de su carbono"
    return None


def _revisar_metoxibenceno() -> str | None:
    """El -O-CH3 del metoxibenceno sale por la derecha del anillo."""
    mol, _, _, _ = _colocada("COc1ccccc1")
    conf = mol.GetConformer()
    o = _unico(mol, 8)
    anillo = postura.anillos_aromaticos(mol)[0]
    borde = max(conf.GetAtomPosition(i).x for i in anillo)
    if conf.GetAtomPosition(o).x <= borde + 0.5:
        return "el O del metoxibenceno no sale por la derecha del anillo"
    return None


def _revisar_ester_diarilico() -> str | None:
    """En el benzoato de fenilo, un anillo a cada lado del rótulo."""
    mol, _, anclas, _ = _colocada("O=C(Oc1ccccc1)c1ccccc1", con_grupos=True)
    conf = mol.GetConformer()
    centro = [i for i, ancla in anclas.items() if ancla == "centro"][0]
    x_centro = conf.GetAtomPosition(centro).x
    centros = [
        sum(conf.GetAtomPosition(i).x for i in anillo) / len(anillo)
        for anillo in postura.anillos_aromaticos(mol)
    ]
    if len(centros) != 2 or (centros[0] - x_centro) * (centros[1] - x_centro) >= 0:
        return "los dos anillos del benzoato de fenilo no quedan a los lados"
    return None


GEOMETRIAS = [
    ("pentenol: el O debajo", _revisar_pentenol),
    ("metoxibenceno: el O a la derecha", _revisar_metoxibenceno),
    ("éster diarílico: un anillo a cada lado", _revisar_ester_diarilico),
]


def main() -> int:
    fallos: list[str] = []

    for nombre, smiles, con_grupos, esperada in CASOS:
        _, _, _, registro = _colocada(smiles, con_grupos)
        if registro.regla != esperada:
            fallos.append(
                f"{nombre}: se esperaba la regla «{esperada}», ha salido {registro.regla!r}"
            )

    for nombre, revisar in GEOMETRIAS:
        problema = revisar()
        if problema:
            fallos.append(f"{nombre}: {problema}")

    print(f"Casos de regla: {len(CASOS)}   geometrías: {len(GEOMETRIAS)}")
    print("Reglas de postura comprobadas:")
    for regla in [r.nombre for r in postura.REGLAS]:
        casos = sum(1 for _, _, _, esperada in CASOS if esperada == regla)
        print(f"  · {regla} ({casos})")

    if fallos:
        print(f"\nFALLOS ({len(fallos)}):")
        for f in fallos:
            print("  -", f)
        return 1

    print("\nCada molécula se coloca con su regla, donde dice la postura.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
