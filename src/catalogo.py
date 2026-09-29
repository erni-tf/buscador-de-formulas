"""Carga y validación de la base de datos de compuestos.

La base vive en `datos/compuestos.toml` y se edita a mano.
Este módulo la lee, la valida y la deja lista para el resto del programa.
"""

from __future__ import annotations

import tomllib
from dataclasses import dataclass, field
from pathlib import Path

RAIZ = Path(__file__).resolve().parent.parent
RUTA_DATOS = RAIZ / "datos" / "compuestos.toml"

# Familias válidas (para el filtro de la aplicación)
FAMILIAS = [
    "alcanos",
    "alquenos",
    "alquinos",
    "aromáticos",
    "alcoholes",
    "éteres",
    "cetonas",
    "aldehídos",
    "ácidos",
    "aminas",
    "ésteres",
    "amidas",
    "sustituyentes",
    "grupos funcionales",
    "conceptos",
]


@dataclass
class Entrada:
    """Un compuesto o una ficha de la base de datos."""

    id: str
    nombre: str
    tipo: str  # "compuesto" o "ficha"
    alias: list[str] = field(default_factory=list)
    familias: list[str] = field(default_factory=list)
    smiles: str = ""
    condensada: str = ""
    notas: str = ""

    @property
    def familia(self) -> str:
        """Familia principal (la primera)."""
        return self.familias[0] if self.familias else ""

    @property
    def tiene_anillo(self) -> bool:
        """True si la molécula tiene algún ciclo (se dibuja con RDKit)."""
        from rdkit import Chem

        mol = Chem.MolFromSmiles(self.smiles)
        if mol is None:
            return True
        return mol.GetRingInfo().NumRings() > 0

    @property
    def todos_los_nombres(self) -> list[str]:
        return [self.nombre] + list(self.alias)


def cargar(ruta: Path = RUTA_DATOS) -> list[Entrada]:
    """Lee el fichero TOML y devuelve la lista de entradas."""
    with open(ruta, "rb") as f:
        datos = tomllib.load(f)

    entradas: list[Entrada] = []
    for tipo, clave in (("compuesto", "compuesto"), ("ficha", "ficha")):
        for bruto in datos.get(clave, []):
            familias = bruto.get("familias")
            if familias is None:
                familias = [bruto.get("familia", "")]
            entradas.append(
                Entrada(
                    id=bruto["id"],
                    nombre=bruto["nombre"],
                    tipo=tipo,
                    alias=list(bruto.get("alias", [])),
                    familias=list(familias),
                    smiles=bruto.get("smiles", ""),
                    condensada=bruto.get("condensada", ""),
                    notas=bruto.get("notas", ""),
                )
            )
    return entradas


def validar(entradas: list[Entrada]) -> list[str]:
    """Comprueba la base y devuelve la lista de problemas encontrados."""
    from rdkit import Chem

    problemas: list[str] = []
    ids: set[str] = set()

    for e in entradas:
        if e.id in ids:
            problemas.append(f"id repetido: {e.id}")
        ids.add(e.id)

        if e.familias:
            for familia in e.familias:
                if familia not in FAMILIAS:
                    problemas.append(f"{e.id}: familia desconocida «{familia}»")
        else:
            problemas.append(f"{e.id}: sin familia")

        if not e.smiles:
            problemas.append(f"{e.id}: falta el SMILES")
        else:
            mol = Chem.MolFromSmiles(e.smiles)
            if mol is None:
                problemas.append(f"{e.id}: SMILES inválido «{e.smiles}»")

        if e.tipo == "compuesto" and not e.condensada and not e.tiene_anillo:
            problemas.append(f"{e.id}: sin fórmula condensada y sin anillo (¿falta algo?)")

    return problemas


if __name__ == "__main__":
    lista = cargar()
    fallos = validar(lista)
    print(f"Entradas cargadas: {len(lista)}")
    compuestos = [e for e in lista if e.tipo == "compuesto"]
    fichas = [e for e in lista if e.tipo == "ficha"]
    print(f"  compuestos: {len(compuestos)}")
    print(f"  fichas:     {len(fichas)}")
    if fallos:
        print(f"\nPROBLEMAS ({len(fallos)}):")
        for p in fallos:
            print("  -", p)
        raise SystemExit(1)
    print("\nBase de datos válida.")
