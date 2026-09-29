"""Postura: cómo se coloca una molécula antes de dibujarla.

La postura es lo que hace que un dibujo quede «como en los apuntes»: la
molécula se gira, se dobla o se refleja antes de que RDKit la pinte.

Tres piezas:

  * las **reglas de postura**: una tabla ordenada de formas nombradas de
    colocar la molécula (el 2-penten-3-ol, el metoxibenceno, el éster
    diarílico…). Cada regla sabe a qué moléculas se aplica (`aplica`, que
    sólo mira y devuelve un plan) y cómo las coloca (`aplicar`, que recibe
    el plan y mueve los átomos). La primera regla que reclama la molécula es
    la que actúa;
  * el **quiebro** de las cadenas que CoordGen ha dibujado rectas (los
    alquinos);
  * el **espejo** que deja el vinilo escrito a la derecha de su unión.

`postura()` aplica las tres cosas y devuelve el registro de lo que hizo, que
es lo que miran los tests. Si ninguna regla reclama la molécula, el dibujo se
queda con la postura por defecto de RDKit.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable

from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Geometry import Point3D

# CoordGen produce esquemas más rectos y regulares (cadenas horizontales),
# como los dibujos a mano de los apuntes: es la geometría de partida que
# después ajustan las reglas.
rdDepictor.SetPreferCoordGen(True)

# Ángulo del quiebro en las cadenas que CoordGen ha dibujado rectas.
_DOBLADO = 60.0

# Caracteres que RDKit no escribe a la manera del curso: los rótulos que los
# llevan se dibujan con el motor de las fórmulas semidesarrolladas, y su
# postura se decide de otra forma (van al lado de su unión).
A_MANO = ("≡", "⁻", "-")


# ---------------------------------------------------------------------------
#  Geometría
# ---------------------------------------------------------------------------

def _girar(
    mol: Chem.Mol,
    centro: tuple[float, float],
    angulo: float,
    indices: list[int] | None = None,
) -> None:
    """Gira las coordenadas (en el plano) alrededor de un punto."""
    conf = mol.GetConformer()
    coseno, seno = math.cos(angulo), math.sin(angulo)
    if indices is None:
        indices = list(range(mol.GetNumAtoms()))
    for i in indices:
        p = conf.GetAtomPosition(i)
        dx, dy = p.x - centro[0], p.y - centro[1]
        conf.SetAtomPosition(
            i,
            Point3D(
                centro[0] + coseno * dx - seno * dy,
                centro[1] + seno * dx + coseno * dy,
                p.z,
            ),
        )


def _voltear_vertical(mol: Chem.Mol) -> None:
    """Refleja el dibujo de arriba abajo."""
    conf = mol.GetConformer()
    for i in range(mol.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, Point3D(p.x, -p.y, p.z))


def _angulo_en(mol: Chem.Mol, i: int, j: int, k: int) -> float:
    """Ángulo (en grados) entre los enlaces j-i y j-k."""
    conf = mol.GetConformer()
    p, q, r = (conf.GetAtomPosition(x) for x in (i, j, k))
    ux, uy = p.x - q.x, p.y - q.y
    vx, vy = r.x - q.x, r.y - q.y
    norma = math.hypot(ux, uy) * math.hypot(vx, vy)
    if not norma:
        return 0.0
    coseno = max(-1.0, min(1.0, (ux * vx + uy * vy) / norma))
    return math.degrees(math.acos(coseno))


def _subarbol(mol: Chem.Mol, raiz: int, bloqueado: int) -> list[int]:
    """Átomos que cuelgan de `raiz` sin pasar por `bloqueado`."""
    dentro: list[int] = []
    vistos = {bloqueado}
    pendientes = [raiz]
    while pendientes:
        i = pendientes.pop()
        if i in vistos:
            continue
        vistos.add(i)
        dentro.append(i)
        pendientes.extend(v.GetIdx() for v in mol.GetAtomWithIdx(i).GetNeighbors())
    return dentro


def _direccion(conf: Chem.Conformer, i: int, j: int) -> float:
    """Dirección (ángulo) del enlace i->j."""
    pi, pj = conf.GetAtomPosition(i), conf.GetAtomPosition(j)
    return math.atan2(pj.y - pi.y, pj.x - pi.x)


def _cadenas_rectas(
    mol: Chem.Mol, excluidos: set[int] | None = None
) -> list[list[int]]:
    """Caminos que cuelgan de un anillo y han quedado dibujados casi rectos.

    CoordGen alinea los tramos con carbonos sp (los alquinos) y en el dibujo
    parecen una sola raya larga; en los apuntes cada enlace cambia de
    dirección. Se devuelve el camino completo, empezando por el átomo del
    anillo que lo sujeta, para que el quiebro separe también el primer enlace
    del sustituyente. Los átomos con rótulo propio (los grupos ya escritos,
    como el C(=O)-O del benzoato de fenilo) no cuentan como cadena.
    """
    en_anillo = {i for anillo in mol.GetRingInfo().AtomRings() for i in anillo}
    excluidos = excluidos or set()

    def es_recto(i: int) -> bool:
        atomo = mol.GetAtomWithIdx(i)
        if (
            i in excluidos
            or i in en_anillo
            or len(atomo.GetNeighbors()) != 2
        ):
            return False
        a, b = (v.GetIdx() for v in atomo.GetNeighbors())
        return _angulo_en(mol, a, i, b) >= 170.0

    def vecinos_rectos(i: int) -> list[int]:
        return [
            v.GetIdx()
            for v in mol.GetAtomWithIdx(i).GetNeighbors()
            if es_recto(v.GetIdx())
        ]

    def fuera(i: int, dentro: set[int]) -> list[int]:
        return [
            v.GetIdx()
            for v in mol.GetAtomWithIdx(i).GetNeighbors()
            if v.GetIdx() not in dentro
        ]

    rectos = {a.GetIdx() for a in mol.GetAtoms() if es_recto(a.GetIdx())}
    cadenas: list[list[int]] = []
    vistos: set[int] = set()
    for inicio in sorted(rectos):
        if inicio in vistos:
            continue
        # ir hasta un extremo del tramo de átomos rectos
        extremo, anterior = inicio, None
        while True:
            siguientes = [j for j in vecinos_rectos(extremo) if j != anterior]
            if not siguientes:
                break
            anterior, extremo = extremo, siguientes[0]
        tramo = [extremo]
        vistos.add(extremo)
        while True:
            siguientes = [j for j in vecinos_rectos(tramo[-1]) if j not in vistos]
            if not siguientes:
                break
            tramo.append(siguientes[0])
            vistos.add(siguientes[0])

        dentro = set(tramo)
        if len(tramo) == 1:
            extremos = fuera(tramo[0], dentro)
            if len(extremos) != 2:
                continue
            izquierda, derecha = extremos
        else:
            externos_izq = fuera(tramo[0], dentro)
            externos_der = fuera(tramo[-1], dentro)
            if len(externos_izq) != 1 or len(externos_der) != 1:
                continue
            izquierda, derecha = externos_izq[0], externos_der[0]

        camino = [izquierda] + tramo + [derecha]
        # el camino debe empezar por el lado que cuelga del anillo

        def anillo_de(i: int) -> int | None:
            if i in en_anillo:
                return i
            vecinos = [
                v.GetIdx()
                for v in mol.GetAtomWithIdx(i).GetNeighbors()
                if v.GetIdx() in en_anillo
            ]
            return vecinos[0] if len(vecinos) == 1 else None

        inicio_anillo = anillo_de(camino[0])
        if inicio_anillo is None:
            inicio_anillo = anillo_de(camino[-1])
            if inicio_anillo is None:
                continue
            camino.reverse()
        if camino[0] not in en_anillo:
            camino = [inicio_anillo] + camino
        cadenas.append(camino)
    return cadenas


def _quiebro(mol: Chem.Mol, excluidos: set[int] | None = None) -> int:
    """Da ángulos a las cadenas que han salido rectas (los alquinos).

    Se toma como referencia el enlace que sujeta la cadena al anillo y cada
    enlace siguiente queda a un lado y otro, en zigzag, como en los apuntes.
    Devuelve cuántas cadenas ha doblado.
    """
    dobladas = 0
    for camino in _cadenas_rectas(mol, excluidos):
        conf = mol.GetConformer()
        base = _direccion(conf, camino[0], camino[1])
        for k in range(1, len(camino) - 1):
            conf = mol.GetConformer()
            objetivo = base + math.radians(_DOBLADO if k % 2 else -_DOBLADO)
            actual = _direccion(conf, camino[k], camino[k + 1])
            p = conf.GetAtomPosition(camino[k])
            _girar(
                mol,
                (p.x, p.y),
                objetivo - actual,
                _subarbol(mol, camino[k + 1], camino[k]),
            )
        dobladas += 1
    return dobladas


# ---------------------------------------------------------------------------
#  Anillos aromáticos
# ---------------------------------------------------------------------------

def anillos_aromaticos(mol: Chem.Mol) -> list[tuple[int, ...]]:
    """Anillos en los que todos los átomos son aromáticos (bencenos...)."""
    return [
        anillo
        for anillo in mol.GetRingInfo().AtomRings()
        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in anillo)
    ]


def _anisol(mol: Chem.Mol) -> tuple[tuple[int, ...], int] | None:
    """(anillo, oxígeno) si es un benceno con un único sustituyente -O-CH3.

    Es el éter no simétrico de los apuntes, que se dibuja con el O y el CH3
    separados y con su ralla entre los dos, no con la etiqueta «OCH3».
    """
    anillos = anillos_aromaticos(mol)
    if len(anillos) != 1 or len(anillos[0]) != 6:
        return None
    anillo = anillos[0]
    exociclicos = [
        vecino.GetIdx()
        for i in anillo
        for vecino in mol.GetAtomWithIdx(i).GetNeighbors()
        if vecino.GetIdx() not in anillo
    ]
    if len(exociclicos) != 1:
        return None
    oxigeno = mol.GetAtomWithIdx(exociclicos[0])
    if oxigeno.GetAtomicNum() != 8:
        return None
    tiene_metilo = any(
        v.GetAtomicNum() == 6 and v.GetTotalNumHs() == 3
        for v in oxigeno.GetNeighbors()
        if v.GetIdx() not in anillo
    )
    return (anillo, exociclicos[0]) if tiene_metilo else None


def es_anisol(mol: Chem.Mol) -> bool:
    """True para el metoxibenceno."""
    return _anisol(mol) is not None


# ---------------------------------------------------------------------------
#  Espejo del vinilo y sentido de las etiquetas
# ---------------------------------------------------------------------------

def _vinilos(mol: Chem.Mol, etiquetas: dict[int, str]) -> list[tuple[int, int]]:
    """Pares (átomo del vinilo, átomo al que se une), colapsado o sin colapsar."""
    pares: list[tuple[int, int]] = []
    for i, et in etiquetas.items():
        if et.startswith("CH="):
            vecinos = [v.GetIdx() for v in mol.GetAtomWithIdx(i).GetNeighbors()]
            if vecinos:
                pares.append((i, vecinos[0]))
    if pares:
        return pares
    patron = Chem.MolFromSmarts("[R][CX3H1]=[CX3H2]")
    if patron is None:
        return []
    return [(m[1], m[0]) for m in mol.GetSubstructMatches(patron)]


def _espejo(mol: Chem.Mol, etiquetas: dict[int, str]) -> bool:
    """Refleja el dibujo si el vinilo queda a la izquierda de su unión.

    El grupo se escribe CH=CH2 (con la unión por la izquierda, como en los
    apuntes); si el esqueleto lo coloca a la izquierda del anillo, se refleja
    para que la línea llegue al CH y el texto se lea en el sentido correcto.
    Devuelve True si ha reflejado.
    """
    pares = _vinilos(mol, etiquetas)
    if not pares:
        return False
    conf = mol.GetConformer()
    for i, union in pares:
        p = conf.GetAtomPosition(i)
        q = conf.GetAtomPosition(union)
        if p.x > q.x:
            return False  # el vinilo ya está a la derecha: no hay que hacer nada
    for i in range(mol.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, Point3D(-p.x, p.y, p.z))
    return True


def orientar_etiqueta(mol: Chem.Mol, indice: int, etiqueta: str) -> str:
    """Escribe el grupo en el sentido en que se une a la cadena.

    El etilo es CH2CH3 o CH3CH2 según de qué lado esté el anillo, para que
    se lea de izquierda a derecha en el orden en que se une, como en los
    apuntes (CH3-CH2-anillo).
    """
    if etiqueta != "CH2CH3":
        return etiqueta
    conf = mol.GetConformer()
    p = conf.GetAtomPosition(indice)
    for vecino in mol.GetAtomWithIdx(indice).GetNeighbors():
        q = conf.GetAtomPosition(vecino.GetIdx())
        if q.x > p.x + 0.1:
            return "CH3CH2"  # la cadena está a la derecha
        if q.x < p.x - 0.1:
            return "CH2CH3"  # la cadena está a la izquierda
    return etiqueta


# ---------------------------------------------------------------------------
#  Sustituyentes que van a la derecha
# ---------------------------------------------------------------------------

# Sustituyentes que los apuntes dibujan hacia la derecha; en la vista
# esquelética se reconocen por su SMARTS y en la semidesarrollada por la
# etiqueta (el grupo ya está colapsado).
_SUSTITUYENTES = [
    ("[c][CX2H0]#[CX2H1]", 1),              # etinilbenceno: -C≡CH
    ("[c][CH2][CX2H0]#[CX2H0][CH3]", 1),    # 1-fenil-2-butino: -CH2-C≡C-CH3
    ("[c][CH2][CH2][OX2H1]", 1),            # 2-feniletanol: -CH2-CH2-OH
    ("[c][CH2][NX3H2]", 1),                 # bencilamina: -CH2NH2
    ("[c][CH2][#0]", 1),                    # bencilo: -CH2⁻ (punto de unión)
]


def _atomo_sustituyente(
    mol: Chem.Mol, etiquetas: dict[int, str], anclas: dict[int, int | str] | None = None
) -> int | None:
    """Átomo del sustituyente que hay que colocar a la derecha."""
    anclas = anclas or {}
    for i, etiqueta in etiquetas.items():
        if any(c in etiqueta for c in A_MANO):
            # si el rótulo se une por su último grupo, va hacia la izquierda:
            # no se toca la orientación. Si va entre dos trozos (centro),
            # tampoco.
            return None if anclas.get(i) in ("ultimo", "centro") else i
    for smarts, indice in _SUSTITUYENTES:
        patron = Chem.MolFromSmarts(smarts)
        if patron is None:
            continue
        coincidencias = mol.GetSubstructMatches(patron)
        if coincidencias:
            return coincidencias[0][indice]
    return None


# ---------------------------------------------------------------------------
#  Reglas de postura
# ---------------------------------------------------------------------------

def _aplica_penten_3_ol(
    mol: Chem.Mol, etiquetas: dict[int, str], anclas: dict[int, int | str]
) -> tuple[int, int, int, int, int, int] | None:
    """2-penten-3-ol: la cadena horizontal y el -OH debajo, como los apuntes.

    CoordGen lo dibuja en forma de «V» (los dos metilos quedan juntos); aquí
    se apunta la cadena para recolocarla en zigzag horizontal, de izquierda a
    derecha.
    """
    if Chem.MolToSmiles(mol) != "CC=C(O)CC":
        return None

    def hidrogenos(i: int) -> int:
        return mol.GetAtomWithIdx(i).GetTotalNumHs()

    # C1 (metilo del etilo) -> C2 (CH2) -> C3 (C-OH) = C4 (CH) -> C5 (metilo)
    c1 = c2 = c3 = c4 = c5 = None
    for m in (
        a.GetIdx()
        for a in mol.GetAtoms()
        if a.GetAtomicNum() == 6 and hidrogenos(a.GetIdx()) == 3
    ):
        vecino = mol.GetAtomWithIdx(m).GetNeighbors()[0].GetIdx()
        if hidrogenos(vecino) == 2:      # -CH2-: este metilo es el del etilo
            c1, c2 = m, vecino
            break
    if c1 is None:
        return None
    for v in mol.GetAtomWithIdx(c2).GetNeighbors():
        if v.GetIdx() != c1:
            c3 = v.GetIdx()
    if c3 is None or hidrogenos(c3) != 0:
        return None
    oxigeno = None
    for v in mol.GetAtomWithIdx(c3).GetNeighbors():
        if v.GetAtomicNum() == 8:
            oxigeno = v.GetIdx()
        elif mol.GetBondBetweenAtoms(c3, v.GetIdx()).GetBondType() == Chem.BondType.DOUBLE:
            c4 = v.GetIdx()
    if oxigeno is None or c4 is None:
        return None
    for v in mol.GetAtomWithIdx(c4).GetNeighbors():
        if v.GetIdx() != c3 and hidrogenos(v.GetIdx()) == 3:
            c5 = v.GetIdx()
    if c5 is None:
        return None
    return c1, c2, c3, c4, c5, oxigeno


def _aplicar_penten_3_ol(mol: Chem.Mol, plan: tuple[int, ...]) -> None:
    c1, c2, c3, c4, c5, oxigeno = plan
    conf = mol.GetConformer()
    largo = 1.5
    seno, coseno = math.sin(math.radians(30.0)), math.cos(math.radians(30.0))
    posiciones = {c1: (0.0, 0.0)}
    px, py = posiciones[c1]
    for carbono, arriba in ((c2, True), (c3, False), (c4, True), (c5, False)):
        py = py + (largo * seno if arriba else -largo * seno)
        px = px + largo * coseno
        posiciones[carbono] = (px, py)
    posiciones[oxigeno] = (posiciones[c3][0], posiciones[c3][1] - largo)
    for i, (x, y) in posiciones.items():
        conf.SetAtomPosition(i, Point3D(x, y, 0.0))


def _aplica_metoxibenceno(
    mol: Chem.Mol, etiquetas: dict[int, str], anclas: dict[int, int | str]
) -> tuple[tuple[int, ...], int] | None:
    """Anisol (el éter no simétrico): el -OCH3 a la derecha del anillo.

    Es como está en los apuntes: el anillo a la izquierda y «O-CH3» al lado.
    """
    return _anisol(mol)


def _aplicar_metoxibenceno(mol: Chem.Mol, plan: tuple[tuple[int, ...], int]) -> None:
    anillo, oxigeno = plan
    conf = mol.GetConformer()
    centro = (sum(conf.GetAtomPosition(i).x for i in anillo) / 6)
    centro_y = (sum(conf.GetAtomPosition(i).y for i in anillo) / 6)
    po = conf.GetAtomPosition(oxigeno)
    _girar(mol, (centro, centro_y), -math.atan2(po.y - centro_y, po.x - centro))


def _aplica_ciclobutano(
    mol: Chem.Mol, etiquetas: dict[int, str], anclas: dict[int, int | str]
) -> tuple | None:
    """1-etil-2-metilciclobutano: el anillo cuadrado y los grupos a la derecha.

    Es como lo dibuja la imagen; CoordGen lo deja girado y con el etilo a la
    izquierda. Vale igual para la vista semidesarrollada, donde el etilo ya
    está escrito como etiqueta.
    """
    anillos = mol.GetRingInfo().AtomRings()
    if len(anillos) != 1 or len(anillos[0]) != 4:
        return None
    if any(a.GetAtomicNum() != 6 for a in mol.GetAtoms()):
        return None  # solo hidrocarburos
    anillo = anillos[0]
    sustituidos = [
        i
        for i in anillo
        if any(v.GetIdx() not in anillo for v in mol.GetAtomWithIdx(i).GetNeighbors())
    ]
    if len(sustituidos) != 2:
        return None
    vecinos = {
        frozenset((i, v.GetIdx()))
        for i in anillo
        for v in mol.GetAtomWithIdx(i).GetNeighbors()
        if v.GetIdx() in anillo
    }
    if frozenset(sustituidos) not in vecinos:  # los dos grupos, en carbonos seguidos
        return None
    conf = mol.GetConformer()
    centro = (
        sum(conf.GetAtomPosition(i).x for i in anillo) / 4,
        sum(conf.GetAtomPosition(i).y for i in anillo) / 4,
    )
    medio = (
        sum(conf.GetAtomPosition(i).x for i in sustituidos) / 2,
        sum(conf.GetAtomPosition(i).y for i in sustituidos) / 2,
    )
    # el etilo, reconocido por su etiqueta o por sus dos hidrógenos, y los
    # pares (átomo del anillo, átomo del grupo) que salen hacia fuera
    etilo = None
    pares: list[tuple[int, int]] = []
    for i in sustituidos:
        for v in mol.GetAtomWithIdx(i).GetNeighbors():
            j = v.GetIdx()
            if j in anillo:
                continue
            pares.append((i, j))
            etiqueta = etiquetas.get(j, "")
            if etiqueta.startswith(("CH2CH3", "CH3CH2")) or v.GetTotalNumHs() == 2:
                etilo = i
    return centro, medio, etilo, pares


def _aplicar_ciclobutano(mol: Chem.Mol, plan: tuple) -> None:
    centro, medio, etilo, pares = plan
    _girar(
        mol,
        centro,
        -math.atan2(medio[1] - centro[1], medio[0] - centro[0]),
    )
    # el etilo queda abajo, como en la imagen (el metilo arriba)
    if etilo is not None and mol.GetConformer().GetAtomPosition(etilo).y > centro[1]:
        _voltear_vertical(mol)
    # los enlaces de los dos grupos salen horizontales hacia la derecha, como
    # en la imagen (si no, quedan en diagonal)
    conf = mol.GetConformer()
    for i, j in pares:
        p, q = conf.GetAtomPosition(i), conf.GetAtomPosition(j)
        _girar(
            mol,
            (p.x, p.y),
            -math.atan2(q.y - p.y, q.x - p.x),
            _subarbol(mol, j, i),
        )


def _aplica_ciclohexanocarbaldehido(
    mol: Chem.Mol, etiquetas: dict[int, str], anclas: dict[int, int | str]
) -> tuple[int, int] | None:
    """Ciclohexanocarbaldehído: la línea del -CHO hacia arriba.

    En los apuntes sale del anillo hacia arriba (no hacia la derecha), tanto
    en la vista de estructura como en la semidesarrollada.
    """
    aldehido = None
    for i, etiqueta in etiquetas.items():
        if etiqueta == "CHO":
            aldehido = i
    if aldehido is None:
        patron = Chem.MolFromSmarts("[CX3H1](=O)[CH1]1CCCCC1")
        if patron is not None:
            coincidencias = mol.GetSubstructMatches(patron)
            if coincidencias:
                aldehido = coincidencias[0][0]
    if aldehido is None:
        return None
    en_anillo = {i for anillo in mol.GetRingInfo().AtomRings() for i in anillo}
    uniones = [
        v.GetIdx()
        for v in mol.GetAtomWithIdx(aldehido).GetNeighbors()
        if v.GetIdx() in en_anillo
    ]
    if len(uniones) != 1:
        return None
    return aldehido, uniones[0]


def _aplicar_ciclohexanocarbaldehido(mol: Chem.Mol, plan: tuple[int, int]) -> None:
    aldehido, union = plan
    conf = mol.GetConformer()
    p = conf.GetAtomPosition(aldehido)
    q = conf.GetAtomPosition(union)
    _girar(
        mol,
        (q.x, q.y),
        math.radians(90.0) - math.atan2(p.y - q.y, p.x - q.x),
    )


def _aplica_acido_ciclohexilpropanoico(
    mol: Chem.Mol, etiquetas: dict[int, str], anclas: dict[int, int | str]
) -> tuple[int, int] | None:
    """Ácido 2-ciclohexilpropanoico: el anillo arriba y la cadena debajo.

    Así la línea baja del anillo a la cadena y el CH3-CH-COOH cuelga justo
    debajo, como en los apuntes.
    """
    atomo = None
    for i, etiqueta in etiquetas.items():
        if etiqueta == "CH3-CH-COOH":
            atomo = i
    if atomo is None:
        patron = Chem.MolFromSmarts("[C][CX4H1]([CH3])C(=O)O")
        if patron is not None:
            coincidencias = mol.GetSubstructMatches(patron)
            if coincidencias:
                atomo = coincidencias[0][1]
    if atomo is None:
        return None
    en_anillo = {i for anillo in mol.GetRingInfo().AtomRings() for i in anillo}
    uniones = [
        v.GetIdx()
        for v in mol.GetAtomWithIdx(atomo).GetNeighbors()
        if v.GetIdx() in en_anillo
    ]
    if len(uniones) != 1:
        return None
    return atomo, uniones[0]


def _aplicar_acido_ciclohexilpropanoico(mol: Chem.Mol, plan: tuple[int, int]) -> None:
    atomo, union = plan
    conf = mol.GetConformer()
    p = conf.GetAtomPosition(atomo)
    q = conf.GetAtomPosition(union)
    # la dirección anillo->cadena, hacia abajo
    _girar(
        mol,
        (q.x, q.y),
        math.radians(-90.0) - math.atan2(p.y - q.y, p.x - q.x),
    )


def _aplica_ester_diarilico(
    mol: Chem.Mol, etiquetas: dict[int, str], anclas: dict[int, int | str]
) -> tuple[int, list[int]] | None:
    """Éster entre dos anillos (benzoato de fenilo): los dos anillos a los lados.

    El rótulo C(=O)-O va en medio y cada anillo se gira alrededor del carbono
    del éster para que su enlace quede horizontal (el vértice del anillo
    apuntando al rótulo), como en los apuntes.
    """
    for i, ancla in anclas.items():
        if ancla != "centro" or i not in etiquetas:
            continue
        vecinos = [v.GetIdx() for v in mol.GetAtomWithIdx(i).GetNeighbors()]
        if len(vecinos) != 2:
            continue
        return i, vecinos
    return None


def _aplicar_ester_diarilico(mol: Chem.Mol, plan: tuple[int, list[int]]) -> None:
    centro, vecinos = plan
    conf = mol.GetConformer()
    p = conf.GetAtomPosition(centro)
    for j in vecinos:
        pj = conf.GetAtomPosition(j)
        dx, dy = pj.x - p.x, pj.y - p.y
        objetivo = 0.0 if dx >= 0 else math.pi
        _girar(
            mol,
            (p.x, p.y),
            objetivo - math.atan2(dy, dx),
            _subarbol(mol, j, centro),
        )


def _aplica_sustituyentes(
    mol: Chem.Mol, etiquetas: dict[int, str], anclas: dict[int, int | str]
) -> tuple[int, int] | None:
    """Pone el sustituyente hacia la derecha, como en los apuntes."""
    atomo = _atomo_sustituyente(mol, etiquetas, anclas)
    if atomo is None:
        return None
    en_anillo = {i for anillo in mol.GetRingInfo().AtomRings() for i in anillo}
    uniones = [
        v.GetIdx()
        for v in mol.GetAtomWithIdx(atomo).GetNeighbors()
        if v.GetIdx() in en_anillo
    ]
    if len(uniones) != 1:
        return None
    return atomo, uniones[0]


def _aplicar_sustituyentes(mol: Chem.Mol, plan: tuple[int, int]) -> None:
    atomo, union = plan
    conf = mol.GetConformer()
    p = conf.GetAtomPosition(atomo)
    q = conf.GetAtomPosition(union)
    _girar(mol, (q.x, q.y), -math.atan2(p.y - q.y, p.x - q.x))


# ---------------------------------------------------------------------------
#  La tabla y el pipeline
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Regla:
    """Una forma nombrada de colocar la molécula.

    `aplica` sólo mira: devuelve lo que hace falta para colocarla (el plan),
    o None si la molécula no es de las suyas. `aplicar` recibe ese plan y
    siempre puede terminar.
    """

    nombre: str
    aplica: Callable[..., object | None]
    aplicar: Callable[..., None]


# El orden es la precedencia: la primera regla que reclama la molécula es la
# que la coloca. `sustituyentes` va al final porque es la única general.
REGLAS: list[Regla] = [
    Regla("penten_3_ol", _aplica_penten_3_ol, _aplicar_penten_3_ol),
    Regla("metoxibenceno", _aplica_metoxibenceno, _aplicar_metoxibenceno),
    Regla("ciclobutano", _aplica_ciclobutano, _aplicar_ciclobutano),
    Regla(
        "ciclohexanocarbaldehido",
        _aplica_ciclohexanocarbaldehido,
        _aplicar_ciclohexanocarbaldehido,
    ),
    Regla(
        "acido_ciclohexilpropanoico",
        _aplica_acido_ciclohexilpropanoico,
        _aplicar_acido_ciclohexilpropanoico,
    ),
    Regla("ester_diarilico", _aplica_ester_diarilico, _aplicar_ester_diarilico),
    Regla("sustituyentes", _aplica_sustituyentes, _aplicar_sustituyentes),
]


@dataclass(frozen=True)
class Postura:
    """Lo que la postura ha hecho con una molécula."""

    regla: str | None
    pasos: tuple[str, ...]


def postura(
    mol: Chem.Mol,
    etiquetas: dict[int, str] | None = None,
    anclas: dict[int, int | str] | None = None,
) -> Postura:
    """Coloca la molécula como en los apuntes y cuenta lo que ha hecho."""
    etiquetas = etiquetas or {}
    anclas = anclas or {}
    regla = None
    for r in REGLAS:
        plan = r.aplica(mol, etiquetas, anclas)
        if plan is not None:
            r.aplicar(mol, plan)
            regla = r.nombre
            break
    pasos: list[str] = []
    if _quiebro(mol, set(etiquetas)):
        pasos.append("quiebro")
    if _espejo(mol, etiquetas):
        pasos.append("espejo")
    return Postura(regla=regla, pasos=tuple(pasos))
