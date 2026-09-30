"""Dibujo de estructuras esqueléticas (estilo «libro de texto»).

Usa RDKit para convertir un SMILES en un SVG limpio, en negro sobre blanco,
sin colores y sin dependencias externas.

Dos variantes:
  * svg_esqueleto()            -> solo líneas (como el ejemplo 2)
  * svg_condensado_anillo()    -> metilos como CH3 y grupos funcionales
                                  agrupados con etiqueta (COOH, CHO, NO2...),
                                  como los escribe el profesor en los apuntes

Los anillos aromáticos se pueden dibujar de las dos formas del curso:
con dobles enlaces alternos (por defecto) o con un círculo en medio
(`circulos=True`).

Algunos compuestos llevan ajustes para quedar como en los apuntes; de la
colocación de la molécula (la postura) se encarga `postura.py`. Los rótulos
que RDKit no sabe escribir (el símbolo ≡, la valencia libre ⁻) se dibujan
con el motor de las fórmulas semidesarrolladas, para que las rayas sean las
mismas.
"""

from __future__ import annotations

import math
import re

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Geometry import Point2D

from dibujo_condensada import (
    ASCENSO,
    DESCENSO,
    GROSOR,
    Y_ENLACE,
    Fuente as FuenteTexto,
    Lienzo,
    analizar,
    centros_grupos,
    dibujar,
)
from postura import (
    A_MANO,
    anillos_aromaticos,
    es_anisol,
    orientar_etiqueta,
    postura,
)

ANCHO = 440
ALTO = 320
# El éster entre dos anillos (el benzoato de fenilo) lleva el rótulo
# C(=O)-O entre los dos: con el ancho normal el hueco queda muy justo y los
# enlaces se ven cortos, así que se le da un lienzo más ancho.
ANCHO_ESTER_DIARILICO = 640
GROSOR_LINEA = 2.0  # grosor de las líneas y de los círculos aromáticos

_ETIQUETAS_R = ["R", "R'", "R''"]

# RDKit escribe con su propia fuente y las fórmulas con la del curso; para
# que las letras dibujadas a mano (C≡CH) salgan del mismo tamaño, se corrige
# la diferencia de altura entre las dos fuentes.
_CORRECCION_FUENTE = 0.708 / 0.662
# Separación entre el final del enlace y el texto dibujado a mano.
_HUECO_TEXTO = 0.28

# Fuente del curso (se carga solo si hace falta dibujar texto a mano).
_FUENTE: FuenteTexto | None = None

# Grupos que se agrupan en una sola etiqueta en la vista semidesarrollada.
# Cada entrada es (SMARTS, etiqueta, átomo del match que lleva la etiqueta,
# átomos del match que se borran). El orden importa: los más específicos
# van primero.
_GRUPOS = [
    ("[C][CX4H1]([CH3])C(=O)O", "CH3-CH-COOH", 1, (2, 3, 4, 5), 1),
    ("[CX3](=O)[OX2][OX2H1]", "COOOH", 0, (1, 2, 3)),  # peroxiácido (mCPBA)
    ("[CX3](=O)[OX2H1]", "COOH", 0, (1, 2)),            # ácido carboxílico
    ("[CX3](=O)[O-]", "COO⁻", 0, (1, 2)),               # carboxilato (anión)
    ("[CX3](=O)[NX3H2]", "CONH2", 0, (1, 2)),           # amida
    ("[CX3](=O)Cl", "COCl", 0, (1, 2)),                 # cloruro de ácido
    ("[NX3+](=O)[O-]", "NO2", 0, (1, 2)),               # nitro
    ("[c][CH2][CX3](=O)[OX2][CH3]", "CH3-O-C[=O]-CH2", 1, (2, 3, 4, 5), "ultimo"),
    ("[c][CX3](=O)[OX2][c]", "C[=O]-O", 1, (2, 3), "centro"),  # benzoato de fenilo
    ("[OX2][CH3]", "OCH3", 0, (1,)),                    # metoxi
    ("[R][CH2][CH3]", "CH2CH3", 1, (2,)),               # etilo: como los apuntes
    ("[R][CX3H1]=[CX3H2]", "CH=CH2", 1, (2,)),          # vinilo: como los apuntes
    ("[c][CX2H0]#[CX2H1]", "C≡CH", 1, (2,)),            # etinilbenceno: los apuntes
    ("[c][CH2][CX2H0]#[CX2H0][CH3]", "CH2-C≡C-CH3", 1, (2, 3, 4)),  # 1-fenil-2-butino
    ("[c][CH2][CH2][OX2H1]", "CH2-CH2-OH", 1, (2, 3)),  # 2-feniletanol: los apuntes
    ("[c][CH2][NX3H2]", "CH2NH2", 1, (2,)),              # bencilamina: los apuntes
    ("[c][NX3H1][CH2][CH3]", "NH-CH2-CH3", 1, (2, 3)),   # N-etilanilina: los apuntes
    ("[c][NX3]([CH2][CH3])[CH2][CH3]", "N[CH2[CH3]]-CH2-CH3", 1, (2, 3, 4, 5)),
    ("[c][CH2][#0]", "CH2⁻", 1, (2,)),                  # bencilo: valencia libre
]

# Etiquetas con COOH que se pueden dibujar desarrolladas: es el conmutador
# «COOH / C(=O)-OH» de la aplicación.
_LABELS_COOH = {
    "COOH": "C[=O]-OH",
    "CH3-CH-COOH": "CH3-CH-C[=O]-OH",
}

# Grupos que se agrupan con el conmutador «CH3» en las fórmulas con anillo:
# los metilos (o los etilos) repetidos del mismo átomo, como los escriben los
# apuntes (N(CH3)2 en la ciclohexildimetilamina, C(CH3)3 en el
# tert-butilbenceno...). Cada entrada es (SMARTS, etiqueta, átomo que lleva la
# etiqueta, átomos que se borran). El orden importa: van del más específico
# (tres metilos) al más general.
_GRUPOS_CH3 = [
    ("[CX4H0]([CH3])([CH3])[CH3]", "C(CH3)3", 0, (1, 2, 3)),
    ("[CX4H1]([CH3])[CH3]", "CH(CH3)2", 0, (1, 2)),
    ("[CX4H0]([CH3])[CH3]", "C(CH3)2", 0, (1, 2)),
    ("[NX3H0]([CH3])[CH3]", "N(CH3)2", 0, (1, 2)),
    ("[NX3H0]([CH2][CH3])[CH2][CH3]", "N(CH3CH2)2", 0, (1, 2, 3, 4)),
    ("[OX2]([CH2][CH3])[CH2][CH3]", "O(CH3CH2)2", 0, (1, 2, 3, 4)),
]


def _limpiar_svg(svg: str) -> str:
    """Quita la declaración XML y los espacios de nombres que no se usan."""
    svg = re.sub(r"<\?xml[^>]*\?>", "", svg).strip()
    svg = svg.replace(" xmlns:rdkit='http://www.rdkit.org/xml'", "")
    svg = svg.replace(" xmlns:xlink='http://www.w3.org/1999/xlink'", "")
    svg = svg.replace(" xml:space='preserve'", "")
    svg = svg.replace(" encoding='iso-8859-1'", "")
    return svg


def _hacer_explicito_el_h_del_aldehido(mol: Chem.Mol) -> Chem.Mol:
    """Añade el hidrógeno del grupo aldehído para que se dibuje, como en los
    apuntes (-C(=O)-H), en vez de poner la etiqueta CHO."""
    patron = Chem.MolFromSmarts("[CX3H1]=O")
    if patron is None:
        return mol
    posiciones = [m[0] for m in mol.GetSubstructMatches(patron)]
    if not posiciones:
        return mol
    return Chem.AddHs(mol, onlyOnAtoms=posiciones)


def _colapsar_grupos(
    mol: Chem.Mol, grupos: list[tuple] | None = None
) -> tuple[Chem.Mol, dict[int, str], dict[int, int | str]]:
    """Sustituye cada grupo funcional terminal por un único átomo con etiqueta.

    Devuelve la molécula «de dibujo», el diccionario {índice: etiqueta} y el
    {índice: ancla}, donde el ancla dice qué grupo de la etiqueta se une al
    resto de la molécula (0 = el primero, «ultimo» = el último, «centro» =
    entre dos trozos). `grupos` permite usar una lista distinta de la
    habitual (el metoxibenceno se dibuja con O y CH3 separados, como en los
    apuntes).
    """
    eliminar: dict[int, str | None] = {}
    anclas: dict[int, int | str] = {}
    dueno: dict[int, int] = {}  # átomo borrado -> átomo que lleva la etiqueta
    for grupo in grupos or _GRUPOS:
        smarts, etiqueta, indice, borrar = grupo[:4]
        ancla = grupo[4] if len(grupo) > 4 else 0
        patron = Chem.MolFromSmarts(smarts)
        if patron is None:
            continue
        for match in mol.GetSubstructMatches(patron):
            if any(i in eliminar for i in match):
                continue  # ya colapsado por otro grupo
            eliminar[match[indice]] = etiqueta
            anclas[match[indice]] = ancla
            for i in borrar:
                eliminar[match[i]] = None
                dueno[match[i]] = match[indice]

    rw = Chem.RWMol()
    mapa: dict[int, int] = {}
    for atomo in mol.GetAtoms():
        i = atomo.GetIdx()
        if i in eliminar and eliminar[i] is None:
            continue
        nuevo = Chem.Atom(atomo.GetAtomicNum())
        nuevo.SetFormalCharge(0 if i in eliminar else atomo.GetFormalCharge())
        mapa[i] = rw.AddAtom(nuevo)

    for enlace in mol.GetBonds():
        a, b = enlace.GetBeginAtomIdx(), enlace.GetEndAtomIdx()
        # los enlaces que salían de un átomo borrado pasan al que lleva la
        # etiqueta (el O del éster del benzoato de fenilo une su anillo al C)
        a = dueno.get(a, a)
        b = dueno.get(b, b)
        if a not in mapa or b not in mapa or a == b:
            continue
        if rw.GetBondBetweenAtoms(mapa[a], mapa[b]) is not None:
            continue
        rw.AddBond(mapa[a], mapa[b], enlace.GetBondType())

    dibujo = rw.GetMol()
    Chem.SanitizeMol(dibujo)
    etiquetas = {mapa[i]: et for i, et in eliminar.items() if et}
    mapa_anclas = {mapa[i]: a for i, a in anclas.items()}
    return dibujo, etiquetas, mapa_anclas


def _sin_dobles_en_anillos(mol: Chem.Mol, anillos: list[tuple[int, ...]]) -> Chem.Mol:
    """Copia de dibujo con los enlaces de esos anillos como líneas simples.

    Para la versión con círculo: el anillo se dibuja con líneas simples y el
    círculo se añade luego en el centro.
    """
    rw = Chem.RWMol(mol)
    for anillo in anillos:
        for k in range(len(anillo)):
            enlace = rw.GetBondBetweenAtoms(anillo[k], anillo[(k + 1) % len(anillo)])
            if enlace is not None:
                enlace.SetBondType(Chem.BondType.SINGLE)
                enlace.SetIsAromatic(False)
    return rw.GetMol()


def _svg_circulos(d: rdMolDraw2D.MolDraw2D, mol: Chem.Mol,
                  anillos: list[tuple[int, ...]]) -> str:
    """Círculos aromáticos, centrados en cada anillo, en coordenadas del SVG."""
    conf = mol.GetConformer()
    partes: list[str] = []
    for anillo in anillos:
        puntos = []
        for i in anillo:
            p = conf.GetAtomPosition(i)
            q = d.GetDrawCoords(Point2D(p.x, p.y))
            puntos.append((q.x, q.y))
        cx = sum(p[0] for p in puntos) / len(puntos)
        cy = sum(p[1] for p in puntos) / len(puntos)
        radio = 0.72 * min(math.hypot(x - cx, y - cy) for x, y in puntos)
        partes.append(
            f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{radio:.1f}" '
            f'style="fill:none;stroke:#000000;stroke-width:{GROSOR_LINEA:.1f}px"/>'
        )
    return "".join(partes)


def _con_subindices(etiqueta: str) -> str:
    """Pasa los dígitos de una etiqueta a subíndice: CH3 -> CH<sub>3</sub>.

    RDKit dibuja las etiquetas de los átomos tal cual, así que hay que
    marcarlas para que el 3 del CH3 salga pequeño y bajo, como en los apuntes.
    """
    return re.sub(r"(\d+)", r"<sub>\1</sub>", etiqueta)


def _es_ciclohexanocarbaldehido(mol: Chem.Mol) -> bool:
    """True para el ciclohexano con un solo -CHO.

    En los apuntes («ciclohexano con butaldheído») el grupo se escribe CHO,
    al lado del anillo.
    """
    return Chem.MolToSmiles(mol) == "O=CC1CCCCC1"


def _fuente_del_curso() -> FuenteTexto:
    """Fuente de las fórmulas, compartida entre dibujos."""
    global _FUENTE
    if _FUENTE is None:
        _FUENTE = FuenteTexto()
    return _FUENTE


def _texto_propio(
    d: rdMolDraw2D.MolDraw2D,
    mol: Chem.Mol,
    indice: int,
    texto: str,
    ancla: int | str = 0,
) -> tuple[str, float, float, float, float]:
    """Dibuja una etiqueta con el motor de las fórmulas semidesarrolladas.

    RDKit no sabe escribir el símbolo ≡, el subíndice ni el ⁻; estos grupos
    se dibujan con el mismo motor y las mismas rayas que las fórmulas
    semidesarrolladas, al lado del enlace (y con los dos puntitos de la N si
    lleva nitrógeno). Devuelve el SVG y los bordes izquierdo y derecho.
    """
    fuente = _fuente_del_curso()
    conf = mol.GetConformer()
    p = conf.GetAtomPosition(indice)
    q = d.GetDrawCoords(Point2D(p.x, p.y))
    tamano = d.FontSize() * _CORRECCION_FUENTE

    lienzo = Lienzo()
    principal = texto.replace("⁻", "")
    cadena = analizar(principal.replace("≡", "#")) if principal else []
    if cadena:
        dibujar(lienzo, fuente, cadena, 0.0, 0.0, tamano)
    # valencia libre (CH2⁻): una rayita pequeña y elevada
    borde = lienzo.x_max if lienzo.elementos else 0.0
    for _ in range(texto.count("⁻")):
        x0 = borde + 0.05 * tamano
        x1 = x0 + 0.30 * tamano
        y = -0.45 * tamano
        lienzo.linea(x0, y, x1, y)
        borde = x1
    primer_grupo = re.split(r"[-\[]", texto, 1)[0]
    if "N" in primer_grupo:
        # los dos puntitos del nitrógeno, encima de su letra
        x_n = fuente.ancho(primer_grupo[: primer_grupo.index("N")], tamano)
        x_n += fuente.ancho("N", tamano) / 2
        lienzo.puntos_n(x_n, -tamano * (0.662 + 0.12), tamano)

    # ¿por dónde cuelga el resto de la molécula?
    vecinos = [v.GetIdx() for v in mol.GetAtomWithIdx(indice).GetNeighbors()]
    direccion = "izquierda"
    if len(vecinos) == 1:
        pv = conf.GetAtomPosition(vecinos[0])
        qv = d.GetDrawCoords(Point2D(pv.x, pv.y))
        dx, dy = qv.x - q.x, qv.y - q.y
        if abs(dx) >= abs(dy):
            direccion = "izquierda" if dx < 0 else "derecha"
        else:
            direccion = "arriba" if dy < 0 else "abajo"

    centros = centros_grupos(fuente, cadena, tamano) if cadena else [0.0]
    if ancla == "centro" and cadena:
        # rótulo entre dos trozos (el éster entre los dos anillos del benzoato
        # de fenilo): se centra en el átomo y los enlaces se recortan luego
        tx = q.x - (lienzo.x_min + lienzo.x_max) / 2
        ty = q.y + Y_ENLACE * tamano
    elif ancla == "ultimo" and cadena and direccion == "derecha":
        # el grupo que se une es el último: la etiqueta sale hacia la izquierda
        tx = q.x - tamano * _HUECO_TEXTO - lienzo.x_max
        ty = q.y + Y_ENLACE * tamano
    elif isinstance(ancla, int) and ancla > 0 and cadena:
        # el grupo que se une está en medio: la etiqueta cuelga debajo
        tx = q.x - centros[ancla]
        if direccion == "abajo":
            ty = q.y - tamano * (DESCENSO + 0.25)
        else:
            ty = q.y + tamano * ASCENSO
    else:
        tx = q.x + tamano * _HUECO_TEXTO
        ty = q.y + Y_ENLACE * tamano

    puntos = ""
    if lienzo.puntos:
        puntos = '<g class="puntos-n" fill="#000000">' + "".join(lienzo.puntos) + "</g>"
    elementos = (
        f'<g transform="translate({tx:.1f},{ty:.1f})" stroke="#000000" '
        f'stroke-width="{tamano * GROSOR:.2f}" stroke-linecap="round" '
        f'fill="#000000">' + "".join(lienzo.elementos) + puntos + "</g>"
    )
    arriba = ty + (lienzo.y_min if lienzo.elementos else 0.0)
    abajo = ty + (lienzo.y_max if lienzo.elementos else 0.0)
    return elementos, tx + lienzo.x_min, tx + max(lienzo.x_max, borde), arriba, abajo


def _ensanchar(
    svg: str,
    mol: Chem.Mol,
    d: rdMolDraw2D.MolDraw2D,
    ancho: int,
    alto: int,
    borde_izq: float,
    borde_der: float,
    borde_arr: float,
    borde_aba: float,
) -> str:
    """Amplía el lienzo para que quepa el texto dibujado a mano.

    Si el texto se sale de la caja, se ensancha el dibujo y se centra todo
    (molécula y texto), para que no quede descentrado.
    """
    conf = mol.GetConformer()
    xs: list[float] = [borde_izq, borde_der]
    ys: list[float] = [borde_arr, borde_aba]
    for i in range(mol.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        q = d.GetDrawCoords(Point2D(p.x, p.y))
        xs.append(q.x)
        ys.append(q.y)
    izquierda, derecha = min(xs), max(xs)
    arriba, abajo = min(ys), max(ys)
    margen = 20.0
    # se amplía solo lo que se sale de la caja (por lados o por arriba/abajo)
    nuevo_ancho, borde_x = ancho, 0.0
    if derecha + 6 > ancho or izquierda < 6:
        nuevo_ancho = max(derecha - izquierda + 2 * margen, ancho)
        borde_x = izquierda - (nuevo_ancho - (derecha - izquierda)) / 2
    nuevo_alto, borde_y = alto, 0.0
    if abajo + 6 > alto or arriba < 6:
        nuevo_alto = max(abajo - arriba + 2 * margen, alto)
        borde_y = arriba - (nuevo_alto - (abajo - arriba)) / 2
    svg = svg.replace(f"width='{ancho}px'", f"width='{nuevo_ancho:.0f}px'", 1)
    svg = svg.replace(f"height='{alto}px'", f"height='{nuevo_alto:.0f}px'", 1)
    svg = svg.replace(
        f"viewBox='0 0 {ancho} {alto}'",
        f"viewBox='{borde_x:.1f} {borde_y:.1f} {nuevo_ancho:.0f} {nuevo_alto:.0f}'",
        1,
    )
    svg = svg.replace(f"width='{ancho}.0'", f"width='{nuevo_ancho:.1f}'", 1)
    svg = svg.replace(f"height='{alto}.0'", f"height='{nuevo_alto:.1f}'", 1)
    svg = svg.replace(
        "x='0.0' y='0.0'", f"x='{borde_x:.1f}' y='{borde_y:.1f}'", 1
    )
    return svg


def _recortar_triples(svg: str, mol: Chem.Mol, d: rdMolDraw2D.MolDraw2D) -> str:
    """Acorta las líneas de fuera de los triples enlaces.

    RDKit las dibuja tan largas como el enlace; en los vértices cerrados
    (los alquinos que se han doblado) la línea de fuera se cruza con el
    enlace vecino. Se recortan por los dos extremos, como hace RDKit con
    las líneas de dentro de los dobles enlaces.
    """
    conf = mol.GetConformer()
    for k, enlace in enumerate(mol.GetBonds()):
        if enlace.GetBondType() != Chem.BondType.TRIPLE:
            continue
        i, j = enlace.GetBeginAtomIdx(), enlace.GetEndAtomIdx()
        pi, pj = conf.GetAtomPosition(i), conf.GetAtomPosition(j)
        qi = d.GetDrawCoords(Point2D(pi.x, pi.y))
        qj = d.GetDrawCoords(Point2D(pj.x, pj.y))
        largo = math.hypot(qj.x - qi.x, qj.y - qi.y)
        if not largo:
            continue
        ux, uy = (qj.x - qi.x) / largo, (qj.y - qi.y) / largo

        def recortar(coincidencia: re.Match) -> str:
            x1, y1 = float(coincidencia.group(2)), float(coincidencia.group(3))
            x2, y2 = float(coincidencia.group(5)), float(coincidencia.group(6))
            # distancia del punto medio al eje del enlace: distingue la línea
            # central de las de fuera
            mx = (x1 + x2) / 2 - qi.x
            my = (y1 + y2) / 2 - qi.y
            fuera = abs(mx * -uy + my * ux)
            if fuera < 1.5:
                return coincidencia.group(0)
            recorte = min(fuera, largo * 0.35)
            x1, y1 = x1 + ux * recorte, y1 + uy * recorte
            x2, y2 = x2 - ux * recorte, y2 - uy * recorte
            return (
                coincidencia.group(1)
                + f"{x1:.1f},{y1:.1f}"
                + coincidencia.group(4)
                + f"{x2:.1f},{y2:.1f}"
                + coincidencia.group(7)
            )

        opciones = rf"atom-{i} atom-{j}|atom-{j} atom-{i}"
        patron = re.compile(
            rf"(<path class='bond-{k} (?:{opciones})' d='M )"
            rf"([-\d.]+),([-\d.]+)( L )([-\d.]+),([-\d.]+)(')"
        )
        svg = patron.sub(recortar, svg)
    return svg


def _puntos_n(
    svg: str,
    mol: Chem.Mol,
    d: rdMolDraw2D.MolDraw2D,
    saltar: set[int],
    etiquetas: dict[int, str],
) -> str:
    """Añade los dos puntitos del nitrógeno (el par de electrones).

    Van en un grupo conmutables (class="puntos-n"): la aplicación los enseña
    o los esconde. Solo se dibujan en los nitrógenos con par libre (aminas,
    amidas): ni en el nitro ni en los nitrilos.
    """
    conf = mol.GetConformer()
    circulos: list[str] = []
    for atomo in mol.GetAtoms():
        i = atomo.GetIdx()
        if i in saltar or atomo.GetAtomicNum() != 7 or atomo.GetFormalCharge() != 0:
            continue
        if etiquetas.get(i, "").startswith("NO"):
            continue  # el nitro no lleva puntitos
        if any(e.GetBondType() != Chem.BondType.SINGLE for e in atomo.GetBonds()):
            continue
        p = conf.GetAtomPosition(i)
        q = d.GetDrawCoords(Point2D(p.x, p.y))
        # la letra de la etiqueta más cercana al átomo es su N
        cajas: list[tuple[float, float, float, float]] = []
        for m in re.finditer(rf"<path class='atom-{i}' d='([^']*)'", svg):
            nums = [float(v) for v in re.findall(r"-?\d+\.?\d*", m.group(1))]
            xs, ys = nums[0::2], nums[1::2]
            if xs and ys:
                cajas.append((min(xs), min(ys), max(xs), max(ys)))
        if not cajas:
            continue

        def distancia(caja: tuple[float, float, float, float]) -> float:
            cx = (caja[0] + caja[2]) / 2
            cy = (caja[1] + caja[3]) / 2
            return (cx - q.x) ** 2 + (cy - q.y) ** 2

        # En las etiquetas propias (N(CH3)2, CONH2...) los glifos van en el
        # orden del texto: la N es el que toca. En las de RDKit (NH2...) la
        # letra del símbolo es la que queda más cerca del átomo.
        etiqueta = re.sub(r"<[^>]+>", "", etiquetas.get(i, ""))
        if "N" in etiqueta and len(cajas) == len(etiqueta):
            caja = cajas[etiqueta.index("N")]
        else:
            caja = min(cajas, key=distancia)
        tamano = d.FontSize()
        cx = (caja[0] + caja[2]) / 2
        cy = caja[1] - tamano * 0.10
        radio = tamano * 0.05
        for signo in (-1, 1):
            x = cx + signo * tamano * 0.115
            circulos.append(
                f'<circle cx="{x:.1f}" cy="{cy:.1f}" r="{radio:.1f}" '
                f'style="stroke:none"/>'
            )
    if not circulos:
        return svg
    grupo = '<g class="puntos-n" fill="#000000">' + "".join(circulos) + "</g>"
    return svg.replace("</svg>", grupo + "</svg>")


def _recortar_enlaces_etiqueta(
    svg: str,
    mol: Chem.Mol,
    d: rdMolDraw2D.MolDraw2D,
    indice: int,
    izquierda: float,
    derecha: float,
) -> str:
    """Acorta los enlaces que llegan a un rótulo dibujado entre dos trozos.

    El enlace acaba en el centro del átomo, así que sin recortarlo se mete
    por debajo del texto del rótulo (el C(=O)-O del benzoato de fenilo).
    """
    conf = mol.GetConformer()
    p = conf.GetAtomPosition(indice)
    q = d.GetDrawCoords(Point2D(p.x, p.y))
    hueco = d.FontSize() * _CORRECCION_FUENTE * 0.12
    for vecino in mol.GetAtomWithIdx(indice).GetNeighbors():
        j = vecino.GetIdx()
        enlace = mol.GetBondBetweenAtoms(indice, j)
        if enlace is None:
            continue
        k = enlace.GetIdx()
        pj = conf.GetAtomPosition(j)
        qj = d.GetDrawCoords(Point2D(pj.x, pj.y))
        lado_izquierdo = qj.x < q.x
        opciones = f"atom-{indice} atom-{j}|atom-{j} atom-{indice}"
        patron = re.compile(
            rf"(<path class='bond-{k} (?:{opciones})' d='M )"
            rf"([-\d.]+),([-\d.]+)( L )([-\d.]+),([-\d.]+)(')"
        )

        def recortar(coincidencia: re.Match) -> str:
            x1, y1 = float(coincidencia.group(2)), float(coincidencia.group(3))
            x2, y2 = float(coincidencia.group(5)), float(coincidencia.group(6))
            objetivo = izquierda - hueco if lado_izquierdo else derecha + hueco
            # se recorta el extremo que está junto al rótulo
            d1 = (x1 - q.x) ** 2 + (y1 - q.y) ** 2
            d2 = (x2 - q.x) ** 2 + (y2 - q.y) ** 2
            if d1 <= d2:
                return (
                    coincidencia.group(1)
                    + f"{objetivo:.1f},{y1:.1f}"
                    + coincidencia.group(4)
                    + f"{x2:.1f},{y2:.1f}"
                    + coincidencia.group(7)
                )
            return (
                coincidencia.group(1)
                + f"{x1:.1f},{y1:.1f}"
                + coincidencia.group(4)
                + f"{objetivo:.1f},{y2:.1f}"
                + coincidencia.group(7)
            )

        svg = patron.sub(recortar, svg)
    return svg


def _dibujar(
    mol: Chem.Mol,
    etiquetas: dict[int, str],
    ancho: int,
    alto: int,
    *,
    metilos: bool,
    circulos: bool = False,
    anclas: dict[int, int | str] | None = None,
) -> str:
    anillos: list[tuple[int, ...]] = []
    if circulos:
        anillos = anillos_aromaticos(mol)
        if anillos:
            mol = _sin_dobles_en_anillos(mol, anillos)

    d = rdMolDraw2D.MolDraw2DSVG(ancho, alto)
    o = d.drawOptions()
    o.clearBackground = True
    o.useBWAtomPalette()
    o.bondLineWidth = GROSOR_LINEA
    # La letra se escala con el dibujo (como en los libros): así una etiqueta
    # (Br, OH, CH3...) no se queda diminuta en una molécula grande.
    o.fixedFontSize = -1
    o.minFontSize = 13
    o.maxFontSize = 34
    o.padding = 0.10
    # En la vista de líneas no se escriben CH3/H3C (solo la semidesarrollada
    # lleva esos rótulos, y allí se añaden a mano más abajo).
    o.explicitMethyl = False
    o.dummiesAreAttachments = True

    # El símbolo ≡ y la valencia libre (⁻) no están en la fuente de RDKit:
    # esos grupos (C≡CH, CH2-C≡C-CH3, CH2⁻...) se dibujan a mano con el motor
    # de las fórmulas semidesarrolladas, al lado del enlace, como en los apuntes.
    anclas = anclas or {}
    textos = {
        i: (et, anclas.get(i, 0))
        for i, et in etiquetas.items()
        if any(c in et for c in A_MANO)
    }

    if metilos:
        for atomo in mol.GetAtoms():
            if (
                atomo.GetAtomicNum() == 6
                and atomo.GetDegree() == 1
                and atomo.GetTotalNumHs() == 3
                and atomo.GetIdx() not in textos
            ):
                o.atomLabels[atomo.GetIdx()] = _con_subindices("CH3")

    # Se prepara el dibujo (kekuliza, calcula coordenadas…), se coloca la
    # molécula como en los apuntes (el registro que devuelve la postura lo
    # usan los tests) y se ajustan las etiquetas que dependen de la
    # orientación.
    mol = rdMolDraw2D.PrepareMolForDrawing(mol)
    postura(mol, etiquetas, anclas)
    for i, etiqueta in etiquetas.items():
        if i in textos:
            o.atomLabels[i] = ""  # el rótulo lo dibuja el motor propio
        else:
            o.atomLabels[i] = _con_subindices(orientar_etiqueta(mol, i, etiqueta))

    d.DrawMolecule(mol)
    d.FinishDrawing()

    svg = _limpiar_svg(d.GetDrawingText())
    svg = _recortar_triples(svg, mol, d)
    if textos:
        piezas: list[str] = []
        borde_izq = borde_arr = float("inf")
        borde_der = borde_aba = float("-inf")
        for i, (texto, ancla) in sorted(textos.items()):
            elemento, izquierda, derecha, arriba, abajo = _texto_propio(
                d, mol, i, texto, ancla
            )
            piezas.append(elemento)
            if ancla == "centro":
                # los enlaces no pueden acabar bajo el texto del rótulo
                svg = _recortar_enlaces_etiqueta(svg, mol, d, i, izquierda, derecha)
            borde_izq = min(borde_izq, izquierda)
            borde_der = max(borde_der, derecha)
            borde_arr = min(borde_arr, arriba)
            borde_aba = max(borde_aba, abajo)
        svg = svg.replace("</svg>", "".join(piezas) + "</svg>")
        if (
            borde_der + 6 > ancho
            or borde_izq < 6
            or borde_aba + 6 > alto
            or borde_arr < 6
        ):
            svg = _ensanchar(
                svg, mol, d, ancho, alto, borde_izq, borde_der, borde_arr, borde_aba
            )
    svg = _puntos_n(svg, mol, d, set(textos), etiquetas)
    if anillos:
        svg = svg.replace("</svg>", _svg_circulos(d, mol, anillos) + "</svg>")
    return svg


def _es_ciclobutino(mol: Chem.Mol) -> bool:
    """True para el ciclobutino: un anillo de 4 carbonos con un triple enlace."""
    if mol.GetNumAtoms() != 4 or mol.GetNumBonds() != 4:
        return False
    if any(a.GetAtomicNum() != 6 for a in mol.GetAtoms()):
        return False
    anillos = mol.GetRingInfo().AtomRings()
    if len(anillos) != 1 or len(anillos[0]) != 4:
        return False
    return sum(1 for e in mol.GetBonds() if e.GetBondType() == Chem.BondType.TRIPLE) == 1


def _svg_ciclobutino() -> str:
    """Ciclobutino dibujado a mano: cuadrado sin girar, como en los apuntes.

    RDKit lo dibuja girado (rombo) y las líneas de dentro del triple enlace
    se juntan con los lados del anillo; aquí van dentro y acortadas.
    """
    lado = 240.0
    x0 = (ANCHO - lado) / 2
    y0 = (ALTO - lado) / 2
    x1, y1 = x0 + lado, y0 + lado
    hueco = lado * 0.12    # separación entre las líneas del triple enlace
    recorte = lado * 0.17  # cuánto se acortan las líneas de dentro y de fuera
    lineas = [
        (x0, y0, x1, y0),          # lado de arriba
        (x0, y1, x1, y1),          # lado de abajo
        (x0, y0, x0, y1),          # lado izquierdo
        (x1, y0, x1, y1),          # lado derecho: línea central del triple enlace
        (x1 - hueco, y0 + recorte, x1 - hueco, y1 - recorte),  # línea de dentro
        (x1 + hueco, y0 + recorte, x1 + hueco, y1 - recorte),  # línea de fuera
    ]
    cuerpo = "".join(
        f"<path d='M {a:.1f},{b:.1f} L {c:.1f},{d:.1f}' style='fill:none;"
        f"stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;"
        f"stroke-linejoin:miter;stroke-opacity:1' />"
        for a, b, c, d in lineas
    )
    return (
        f"<svg version='1.1' baseProfile='full' xmlns='http://www.w3.org/2000/svg' "
        f"width='{ANCHO}px' height='{ALTO}px' viewBox='0 0 {ANCHO} {ALTO}'>"
        f"<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='{ANCHO}.0' "
        f"height='{ALTO}.0' x='0.0' y='0.0'> </rect>"
        f"{cuerpo}</svg>"
    )


def svg_esqueleto(
    smiles: str,
    ancho: int = ANCHO,
    alto: int = ALTO,
    *,
    etiquetas_dummy: bool = False,
    circulos: bool = False,
) -> str:
    """SVG de la estructura esquelética (solo líneas y heteroátomos)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"SMILES inválido: {smiles!r}")

    if _es_ciclobutino(mol):
        return _svg_ciclobutino()

    if etiquetas_dummy:
        etiquetas = {}
        k = 0
        for a in mol.GetAtoms():
            if a.GetAtomicNum() == 0:
                etiquetas[a.GetIdx()] = _ETIQUETAS_R[min(k, len(_ETIQUETAS_R) - 1)]
                k += 1
    else:
        etiquetas = {}

    return _dibujar(mol, etiquetas, ancho, alto, metilos=False, circulos=circulos)


def svg_condensado_anillo(
    smiles: str,
    ancho: int = ANCHO,
    alto: int = ALTO,
    *,
    circulos: bool = False,
    cooh_separado: bool = False,
    ch3_agrupados: bool = False,
) -> str:
    """SVG semidesarrollado para compuestos con anillo (CH3, COOH, CHO...).

    Con `cooh_separado` los grupos COOH se escriben desarrollados
    (C(=O)-OH) en vez de juntos y con `ch3_agrupados` los metilos (o etilos)
    repetidos del mismo átomo se juntan en un grupo (N(CH3)2, C(CH3)3...):
    son los conmutadores de la aplicación.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"SMILES inválido: {smiles!r}")
    if _es_ciclobutino(mol):
        return _svg_ciclobutino()
    grupos = None
    if _es_ciclohexanocarbaldehido(mol):
        # el «ciclohexano con butaldheído» acaba en CHO, como en los apuntes
        grupos = _GRUPOS + [("[CX3H1]=O", "CHO", 0, (1,))]
    else:
        mol = _hacer_explicito_el_h_del_aldehido(mol)
    # En el éter no simétrico, los apuntes escriben O y CH3 separados (con su
    # ralla entre los dos), no la etiqueta «OCH3».
    if es_anisol(mol):
        grupos = [g for g in _GRUPOS if g[1] != "OCH3"]
    if ch3_agrupados:
        base = grupos if grupos is not None else _GRUPOS
        # la etiqueta suelta de la N,N-dietilanilina se sustituye por la
        # agrupada: N(CH3CH2)2
        base = [g for g in base if g[1] != "N[CH2[CH3]]-CH2-CH3"]
        grupos = _GRUPOS_CH3 + base
    dibujo, etiquetas, anclas = _colapsar_grupos(mol, grupos)
    if cooh_separado:
        etiquetas = {i: _LABELS_COOH.get(et, et) for i, et in etiquetas.items()}
    if "centro" in anclas.values():
        # el rótulo va entre dos trozos: hace falta más sitio para que no
        # quede apretado (el benzoato de fenilo)
        ancho = max(ancho, ANCHO_ESTER_DIARILICO)
    return _dibujar(
        dibujo, etiquetas, ancho, alto, metilos=True, circulos=circulos, anclas=anclas
    )


if __name__ == "__main__":
    from pathlib import Path

    salida = Path("/tmp/opencode/salida_esqueleto")
    salida.mkdir(parents=True, exist_ok=True)
    pruebas = {
        "esq_ciclopropiloctano": ("CCC(C)C(C1CC1)CCCC", svg_esqueleto),
        "esq_fenol": ("Oc1ccccc1", svg_esqueleto),
        "con_trimetilbenceno": ("Cc1cc(C)cc(C)c1", svg_condensado_anillo),
        "con_acidoformilbenzoico": ("O=Cc1cccc(C(=O)O)c1", svg_condensado_anillo),
        "con_dinitrofenol": ("O=[N+]([O-])c1ccc(O)c([N+](=O)[O-])c1", svg_condensado_anillo),
        "con_mcpba": ("OOC(=O)c1cccc(Cl)c1", svg_condensado_anillo),
        "con_anisol": ("COc1ccccc1", svg_condensado_anillo),
        "con_ciclohexanocarbald": ("O=CC1CCCCC1", svg_condensado_anillo),
        "con_etoximetilciclohexano": ("CCOC1CCCC(C)(C)C1", svg_condensado_anillo),
    }
    for nombre, (smi, fn) in pruebas.items():
        svg = fn(smi)
        (salida / f"{nombre}.svg").write_text(svg, encoding="utf-8")
        print(f"{nombre}: {len(svg)} bytes")
