"""Dibujo de proyecciones de Haworth (los azúcares de los apuntes).

Los azúcares no se dibujan como las demás moléculas: en los apuntes van en
**proyección de Haworth**, el anillo hexagonal con el oxígeno en el vértice
de arriba a la derecha, el CH2OH colgando del de arriba a la izquierda y los
-OH hacia arriba o hacia abajo. Eso es lo que distingue los anómeros: en el
α el -OH del C1 mira hacia abajo y en el β hacia arriba.

Dos dibujos:

  * `svg_haworth()`:   un anillo suelto (un monosacárido).
  * `svg_disacarido()`: los dos anillos unidos por el oxígeno del enlace
    O-glicosídico, con el H de los dos carbonos que se unen, como en los
    apuntes (el enlace α baja y el β sube).

Cada azúcar se describe con sus **restos**: qué cuelga de cada carbono del
anillo (C1…C5) y hacia dónde («arriba», «abajo» o «al lado», que es como los
apuntes dibujan el -OH del extremo reductor, sin decir α ni β).

Las letras y las rayas las pone el motor de las fórmulas semidesarrolladas
(`dibujo_condensada`): la misma letra, los mismos subíndices y el mismo
grosor que en el resto del programa.
"""

from __future__ import annotations

import math

from dibujo_condensada import (
    ALTO_GLIFO,
    FS,
    GROSOR,
    MARGEN,
    Fuente,
    Lienzo,
    analizar,
    dibujar,
    medir,
)

# --- medidas (en múltiplos del tamaño de fuente) ---------------------------
RADIO = 2.6         # del centro del anillo a cada vértice
ENLACE = 1.10       # largo de los enlaces que salen del anillo
HUECO = 0.12        # separación entre el final de un enlace y sus letras
HUECO_O = 0.10      # pequeño blanco entre el O del anillo y sus lados
DESNIVEL = 1.60     # cuánto baja (o sube) el enlace glicosídico
HUECO_ANILLOS = 3.00  # espacio para el oxígeno del enlace entre los anillos

# El anillo, en grados: el C1 a la derecha, el oxígeno arriba a la derecha y
# el C5 (el del CH2OH) arriba a la izquierda, como en los apuntes.
_ANGULOS = {
    "C1": 0.0,
    "O": 60.0,
    "C5": 120.0,
    "C4": 180.0,
    "C3": 240.0,
    "C2": 300.0,
}

# Los enlaces del anillo, en orden.
_ENLACES = [
    ("C1", "O"),
    ("O", "C5"),
    ("C5", "C4"),
    ("C4", "C3"),
    ("C3", "C2"),
    ("C2", "C1"),
]

# Un resto: hacia dónde mira su enlace y qué se escribe al final.
Restos = dict[str, tuple[str, str]]

# Los azúcares de los apuntes del 14 de octubre. El α y el β se distinguen
# en el C1: el -OH hacia abajo (α) o hacia arriba (β).
GLUCOSA_ALFA: Restos = {
    "C1": ("abajo", "OH"),
    "C2": ("abajo", "OH"),
    "C3": ("arriba", "OH"),
    "C4": ("abajo", "OH"),
    "C5": ("arriba", "CH2OH"),
}

GALACTOSA_BETA: Restos = {
    "C1": ("arriba", "OH"),
    "C2": ("abajo", "OH"),
    "C3": ("arriba", "OH"),
    "C4": ("arriba", "OH"),
    "C5": ("arriba", "CH2OH"),
}

# La glucosa del extremo reductor de un disacárido: el -OH de su carbono
# anomérico sale al lado, como lo dibujan los apuntes (es el extremo libre:
# no se dice si es α o β).
GLUCOSA_REDUCTORA: Restos = {**GLUCOSA_ALFA, "C1": ("al lado", "OH")}


# ---------------------------------------------------------------------------
#  Geometría
# ---------------------------------------------------------------------------

def _vertices(centro: tuple[float, float], tamano: float) -> dict[str, tuple[float, float]]:
    """Los seis vértices del anillo, con el oxígeno arriba a la derecha."""
    radio = RADIO * tamano
    return {
        nombre: (
            centro[0] + radio * math.cos(math.radians(grados)),
            centro[1] - radio * math.sin(math.radians(grados)),
        )
        for nombre, grados in _ANGULOS.items()
    }


def _caja_resto(
    fuente: Fuente, texto: str, vertice: tuple[float, float], posicion: str,
    tamano: float,
) -> tuple[float, float, tuple[float, float]]:
    """Dónde van las letras de un resto y dónde acaba su enlace.

    Devuelve el borde izquierdo de las letras, su línea base y el final del
    enlace que las sujeta: hacia arriba o hacia abajo el enlace sale recto del
    vértice, y «al lado» sale horizontal, con las letras centradas en la altura
    del vértice (es el -OH del extremo reductor, que los apuntes escriben al
    lado).
    """
    x, y = vertice
    # El enlace sale del primer átomo del grupo (la C del CH2OH, la O del OH),
    # no del centro de las letras: así la raya no le pasa por encima al
    # subíndice del CH2OH.
    primer_atomo = fuente.ancho(texto[0], tamano) / 2
    largo = ENLACE * tamano
    hueco = HUECO * tamano
    if posicion == "arriba":
        return x - primer_atomo, y - largo - hueco, (x, y - largo)
    if posicion == "abajo":
        base = y + largo + hueco + ALTO_GLIFO * tamano
        return x - primer_atomo, base, (x, y + largo)
    if posicion == "al lado":
        return x + largo + hueco, y + ALTO_GLIFO * tamano / 2, (x + largo, y)
    raise ValueError(f"Un resto no puede mirar «{posicion}»")


HUECO_ROTULO = 0.20  # margen mínimo entre un rótulo apartado y la raya del anillo


def _desplaza_izquierda(
    fuente: Fuente,
    vertices: dict[str, tuple[float, float]],
    carbono: str,
    posicion: str,
    texto: str,
    tamano: float,
) -> float:
    """Cuánto hay que apartar hacia la izquierda un rótulo que choca con la raya.

    El -OH del C4 (el vértice de más a la izquierda de la glucosa y la
    galactosa) cuelga hacia abajo o hacia arriba justo donde pasa la raya del
    anillo (la del C4-C3, si mira hacia abajo; la del C4-C5, si mira hacia
    arriba). La raya es sagrada y no se toca: se aparta el rótulo hacia fuera
    (hacia la izquierda) lo justo para que la raya quede limpia, como en los
    apuntes. El enlace del propio -OH no se mueve: puede que no acabe centrado
    en la `O`, y da igual, en el papel se escribe así.
    """
    x, y = vertices[carbono]
    if posicion not in ("arriba", "abajo"):
        return 0.0
    if abs(x - min(vx for vx, _ in vertices.values())) > 1e-6:
        return 0.0  # solo el vértice de más a la izquierda necesita apartarse
    vecino = "C5" if posicion == "arriba" else "C3"
    if vecino not in vertices:
        return 0.0
    ax, ay = vertices[vecino]
    if abs(ay - y) < 1e-6:
        return 0.0
    izquierda, base, _ = _caja_resto(fuente, texto, vertices[carbono], posicion, tamano)
    y_arriba, y_abajo = base - ALTO_GLIFO * tamano, base
    pendiente = (ax - x) / (ay - y)  # dx/dy de la raya del anillo que baja del vértice
    # Punto de la caja del rótulo donde la raya pasa más a la izquierda.
    y_cerca = y_arriba if pendiente >= 0 else y_abajo
    raya = x + pendiente * (y_cerca - y)
    derecha = izquierda + fuente.ancho(texto, tamano)
    return max(0.0, derecha + HUECO_ROTULO * tamano - raya)


def _lado(posicion: str) -> float:
    """Lado hacia el que mira un resto: hacia abajo, positivo (el eje del SVG)."""
    if posicion == "abajo":
        return 1.0
    if posicion == "arriba":
        return -1.0
    raise ValueError(f"El enlace O-glicosídico no puede mirar «{posicion}»")


def _recortar(
    enlace: tuple[tuple[float, float], tuple[float, float]],
    caja: tuple[float, float, float, float],
) -> tuple[tuple[float, float], tuple[float, float]]:
    """Para un enlace en el borde de la caja de su rótulo.

    El enlace va del primer punto al segundo, que es el centro del rótulo (un
    oxígeno, que se escribe en medio): se devuelve el mismo enlace acortado,
    para que la raya no se meta por debajo de las letras.
    """
    (x1, y1), (x2, y2) = enlace
    dx, dy = x2 - x1, y2 - y1
    pasos = [
        max(x2 - caja[0], caja[2] - x2) / abs(dx) if dx else math.inf,
        max(y2 - caja[1], caja[3] - y2) / abs(dy) if dy else math.inf,
    ]
    s = min(min(pasos), 1.0)
    return (x1, y1), (x2 - s * dx, y2 - s * dy)


def _puente(
    direccion_izquierda: str,
    direccion_derecha: str,
    p_izquierda: tuple[float, float],
    p_derecha: tuple[float, float],
    tamano: float,
) -> tuple[tuple[float, float], tuple[float, float], tuple[float, float]]:
    """El enlace O-glicosídico: por dónde va y dónde se escribe su oxígeno.

    Cada carbono tira del oxígeno hacia su lado (el -OH que ha desaparecido):
    si los dos miran hacia abajo, el oxígeno baja y el enlace va horizontal,
    como el de la maltosa; si miran a lados distintos, el enlace baja (o sube)
    en diagonal y el oxígeno se escribe en medio, como el de la lactosa.
    """
    desnivel = DESNIVEL * tamano
    caida_izquierda = _lado(direccion_izquierda) * desnivel
    caida_derecha = _lado(direccion_derecha) * desnivel
    inicio = (p_izquierda[0], p_izquierda[1] + caida_izquierda)
    fin = (p_derecha[0], p_derecha[1] + caida_derecha)
    medio = ((inicio[0] + fin[0]) / 2, (inicio[1] + fin[1]) / 2)
    return inicio, fin, medio


# ---------------------------------------------------------------------------
#  Dibujo
# ---------------------------------------------------------------------------

def _etiqueta_centrada(
    lienzo: Lienzo, fuente: Fuente, texto: str, x: float, y: float, tamano: float
) -> tuple[float, float, float, float]:
    """Escribe un rótulo centrado en un punto y devuelve su caja.

    Es el oxígeno del anillo y el del enlace O-glicosídico: van en medio de
    las rayas, y con la caja se recortan los enlaces que llegan.
    """
    cadena = analizar(texto)
    ancho = medir(fuente, cadena, tamano).ancho
    alto = ALTO_GLIFO * tamano
    dibujar(lienzo, fuente, cadena, x - ancho / 2, y + alto / 2, tamano)
    hueco = HUECO_O * tamano
    return (
        x - ancho / 2 - hueco,
        y - alto / 2 - hueco,
        x + ancho / 2 + hueco,
        y + alto / 2 + hueco,
    )


def _resto(
    lienzo: Lienzo, fuente: Fuente, vertice: tuple[float, float], posicion: str,
    texto: str, tamano: float, desplaza_x: float = 0.0,
) -> tuple[float, float]:
    """Dibuja un resto (su enlace y sus letras) y devuelve el final del enlace.

    `desplaza_x` aparta las letras hacia la izquierda sin mover el enlace: es
    la salida para los rótulos que chocan con la raya del anillo (ver
    `_desplaza_izquierda`).
    """
    x, base, fin = _caja_resto(fuente, texto, vertice, posicion, tamano)
    dibujar(lienzo, fuente, analizar(texto), x - desplaza_x, base, tamano)
    lienzo.linea(vertice[0], vertice[1], fin[0], fin[1])
    return fin


def _anillo(
    lienzo: Lienzo, fuente: Fuente, centro: tuple[float, float], restos: Restos,
    tamano: float,
) -> dict[str, tuple[float, float]]:
    """Dibuja un anillo con sus restos y devuelve sus vértices."""
    vertices = _vertices(centro, tamano)
    caja = _etiqueta_centrada(lienzo, fuente, "O", *vertices["O"], tamano)
    for a, b in _ENLACES:
        if a == "O":
            a, b = b, a  # el recorte se hace siempre por el lado del oxígeno
        inicio, fin = (vertices[a], vertices[b])
        if b == "O":
            inicio, fin = _recortar((inicio, fin), caja)
        lienzo.linea(inicio[0], inicio[1], fin[0], fin[1])
    for carbono, (posicion, texto) in sorted(restos.items()):
        if carbono not in vertices:
            raise ValueError(f"El anillo no tiene ningún {carbono}")
        _resto(
            lienzo, fuente, vertices[carbono], posicion, texto, tamano,
            _desplaza_izquierda(fuente, vertices, carbono, posicion, texto, tamano),
        )
    return vertices


def svg_haworth(restos: Restos, *, tamano: float = FS) -> str:
    """SVG de un monosacárido en proyección de Haworth."""
    if "O" in restos:
        raise ValueError("El anillo ya lleva su oxígeno: los restos son de C1 a C5")
    fuente = Fuente()
    lienzo = Lienzo()
    _anillo(lienzo, fuente, (0.0, 0.0), restos, tamano)
    return lienzo.svg(grosor=tamano * GROSOR, margen=tamano * MARGEN)


def svg_disacarido(
    izquierda: Restos,
    derecha: Restos,
    *,
    uniones: tuple[str, str] = ("C1", "C4"),
    tamano: float = FS,
) -> str:
    """SVG del disacárido: los dos anillos unidos por su oxígeno.

    `uniones` son el carbono del anillo de la izquierda y el del anillo de la
    derecha que se unen: el enlace O-glicosídico de los apuntes es el 1→4.
    La forma del enlace la decide el azúcar de la izquierda (si su -OH mira
    hacia abajo, el oxígeno baja: es el enlace α; si mira hacia arriba, sube
    y el enlace va en diagonal: es el β). Los dos carbonos que se unen pierden
    su -OH y llevan su H al otro lado, como en los apuntes.
    """
    carbono_izquierda, carbono_derecha = uniones
    for carbono, azucar in (
        (carbono_izquierda, izquierda),
        (carbono_derecha, derecha),
    ):
        if carbono not in azucar:
            raise ValueError(f"Al azúcar le falta el {carbono} del enlace")
    direccion_izquierda = izquierda[carbono_izquierda][0]
    direccion_derecha = derecha[carbono_derecha][0]

    fuente = Fuente()
    lienzo = Lienzo()
    radio = RADIO * tamano
    separacion = 2 * radio + HUECO_ANILLOS * tamano
    sueltos_izquierda = {c: r for c, r in izquierda.items() if c != carbono_izquierda}
    sueltos_derecha = {c: r for c, r in derecha.items() if c != carbono_derecha}
    vertices_izquierda = _anillo(lienzo, fuente, (0.0, 0.0), sueltos_izquierda, tamano)
    vertices_derecha = _anillo(
        lienzo, fuente, (separacion, 0.0), sueltos_derecha, tamano
    )

    # el oxígeno del enlace, en medio de los dos anillos
    inicio, fin, medio = _puente(
        direccion_izquierda,
        direccion_derecha,
        vertices_izquierda[carbono_izquierda],
        vertices_derecha[carbono_derecha],
        tamano,
    )
    caja = _etiqueta_centrada(lienzo, fuente, "O", *medio, tamano)
    for extremo in (inicio, fin):
        (x1, y1), (x2, y2) = _recortar((extremo, medio), caja)
        lienzo.linea(x1, y1, x2, y2)

    # el H de cada carbono del enlace, al otro lado que su -OH, con su raya
    # (la que baja o sube hasta el oxígeno del enlace)
    for vertices, carbono, direccion, extremo in (
        (vertices_izquierda, carbono_izquierda, direccion_izquierda, inicio),
        (vertices_derecha, carbono_derecha, direccion_derecha, fin),
    ):
        contraria = "arriba" if _lado(direccion) > 0 else "abajo"
        final = _resto(lienzo, fuente, vertices[carbono], contraria, "H", tamano)
        lienzo.linea(final[0], final[1], extremo[0], extremo[1])

    return lienzo.svg(grosor=tamano * GROSOR, margen=tamano * MARGEN)


if __name__ == "__main__":
    from pathlib import Path

    salida = Path("/tmp/opencode/salida_haworth")
    salida.mkdir(parents=True, exist_ok=True)
    pruebas = {
        "alfa_glucosa": lambda: svg_haworth(GLUCOSA_ALFA),
        "beta_galactosa": lambda: svg_haworth(GALACTOSA_BETA),
        "maltosa": lambda: svg_disacarido(GLUCOSA_ALFA, GLUCOSA_REDUCTORA),
        "lactosa": lambda: svg_disacarido(GALACTOSA_BETA, GLUCOSA_REDUCTORA),
    }
    for nombre, dibujo in pruebas.items():
        svg = dibujo()
        (salida / f"{nombre}.svg").write_text(svg, encoding="utf-8")
        print(f"{nombre}: {len(svg)} bytes")
