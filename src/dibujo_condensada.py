"""Dibujo de fórmulas semidesarrolladas de cadena abierta (estilo ejemplo 1).

Convierte una fórmula escrita en el DSL propio en un SVG limpio:

    CH3-C{=O}-CH2-CH3

  * grupos:  CH3, CH2, CH, C, O, OH, F, Cl, Br, I, N, NH2, CHO, COOH, NO2...
  * enlaces: - simple, = doble, # triple
  * rama arriba: { ... }      rama abajo: [ ... ]
  * dentro de una rama se puede indicar el enlace:  C{=O}
  * un guion al principio o al final dibuja un enlace abierto (fichas)

El texto se convierte a curvas con fontTools, así el dibujo se ve igual en
cualquier ordenador aunque no tenga la fuente instalada.

La fórmula se puede variar con los conmutadores de la aplicación:
`cooh_separado` dibuja los COOH como C(=O)-OH, `ch2_agrupados` junta los
CH2 seguidos en (CH2)n y `ch3_agrupados` junta en un grupo los CH3 (o los
CH3CH2) que cuelgan del mismo átomo.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

from fontTools.pens.svgPathPen import SVGPathPen
from fontTools.ttLib import TTFont

RUTA_FUENTE = Path("/usr/share/fonts/liberation/LiberationSerif-Regular.ttf")

# --- medidas (en múltiplos del tamaño de fuente) ---------------------------
FS = 46.0            # tamaño de fuente base
FS_SUB = 0.62        # tamaño relativo de los subíndices
SUB_BAJADA = 0.20    # cuánto baja la línea base del subíndice
ENLACE = 0.65        # longitud de un enlace horizontal
# El texto siempre crece hacia arriba desde su línea base, así que una rama
# de arriba y una de abajo no necesitan la misma distancia.
RAMAS_ARRIBA = 1.40  # de la línea base a la línea base de la rama de arriba
RAMAS_ABAJO = 1.40   # de la línea base a la línea base de la rama de abajo
INICIO_ARRIBA = 0.78 # el enlace empieza por encima de las mayúsculas
INICIO_ABAJO = 0.16  # y algo por debajo de la línea base
ALTO_GLIFO = 0.662   # altura de las mayúsculas y los dígitos de la fuente
FIN_RAMA = 0.16      # hueco entre el final del enlace y la rama
HUECO_RAMAS = 0.50   # hueco mínimo entre dos ramas vecinas del mismo lado
GROSOR = 0.075       # grosor de línea
MARGEN = 0.35        # margen alrededor del dibujo
ASCENSO = 0.72       # altura del texto por encima de la línea base
DESCENSO = 0.24      # profundidad por debajo (subíndices y trazos bajos)
Y_ENLACE = 0.22      # altura a la que se dibujan los enlaces simples y dobles
Y_TRIPLE = ALTO_GLIFO / 2  # el triple enlace se centra en la altura de la letra
SEPARACION_DOBLE = 0.07   # separación de las dos líneas del doble enlace
SEPARACION_TRIPLE = 0.22  # separación entre las líneas del triple enlace
ANCHO_SUPER = 0.32        # ancho de la rayita del superíndice (O-)
SUBE_SUPER = 0.45         # cuánto sube el superíndice
RADIO_PUNTO = 0.045       # radio de los dos puntitos del nitrógeno
SEPARACION_PUNTOS = 0.09  # separación entre los dos puntitos


# ---------------------------------------------------------------------------
#  Fuente: medidas y contornos
# ---------------------------------------------------------------------------

class Fuente:
    """Fuente TrueType con la que se miden y dibujan los grupos."""

    def __init__(self, ruta: Path = RUTA_FUENTE) -> None:
        self.tt = TTFont(str(ruta))
        self.gs = self.tt.getGlyphSet()
        self.cmap = self.tt.getBestCmap()
        self.upem = self.tt["head"].unitsPerEm
        self._anchos: dict[tuple[str, float], float] = {}

    def ancho(self, texto: str, tamano: float) -> float:
        """Anchura del texto en píxeles."""
        clave = (texto, round(tamano, 2))
        if clave not in self._anchos:
            total = 0.0
            for caracter in texto:
                glifo = self.cmap.get(ord(caracter))
                if glifo is None:
                    continue
                total += self.gs[glifo].width * tamano / self.upem
            self._anchos[clave] = total
        return self._anchos[clave]

    def contornos(self, texto: str, x: float, y: float, tamano: float) -> list[str]:
        """Contornos SVG de cada carácter, con la línea base en `y`."""
        partes: list[str] = []
        escala = tamano / self.upem
        cursor = x
        for caracter in texto:
            glifo = self.cmap.get(ord(caracter))
            if glifo is None:
                continue
            trazo = self.gs[glifo]
            pluma = SVGPathPen(self.gs, ntos=lambda v: f"{v:.1f}")
            trazo.draw(pluma)
            d = pluma.getCommands()
            if d:
                partes.append(
                    f'<path transform="translate({cursor:.1f},{y:.1f}) '
                    f'scale({escala:.4f},{-escala:.4f})" d="{d}"/>'
                )
            cursor += trazo.width * escala
        return partes


# ---------------------------------------------------------------------------
#  Modelo de la fórmula
# ---------------------------------------------------------------------------

@dataclass
class Trozo:
    """Un trozo de texto de un grupo: letras normales o subíndice."""

    texto: str
    subindice: bool


@dataclass
class Nodo:
    """Un grupo de átomos con sus ramas verticales."""

    grupo: str
    arriba: list["Rama"] = field(default_factory=list)
    abajo: list["Rama"] = field(default_factory=list)

    def trozos(self) -> list[Trozo]:
        """Divide el grupo en texto, subíndices y el signo ⁻: CH3 -> CH + 3."""
        return [
            Trozo(t, t.isdigit())
            for t in re.findall(r"\d+|⁻|[^\d⁻]+", self.grupo)
        ]


@dataclass
class Rama:
    """Una rama vertical: su enlace con el grupo padre y su cadena."""

    enlace: str
    cadena: "Cadena"


Cadena = list[tuple[Nodo, str]]  # (nodo, enlace que le sigue: "-", "=", "#", "")


def analizar(texto: str) -> Cadena:
    """Convierte la fórmula del DSL en una cadena de nodos."""
    pos = 0

    def saltar_espacios() -> None:
        nonlocal pos
        while pos < len(texto) and texto[pos].isspace():
            pos += 1

    def leer_ramas(nodo: Nodo) -> None:
        nonlocal pos
        saltar_espacios()
        while pos < len(texto) and texto[pos] in "{[":
            lado = nodo.arriba if texto[pos] == "{" else nodo.abajo
            cierre = "}" if texto[pos] == "{" else "]"
            pos += 1
            saltar_espacios()
            enlace = "-"
            if pos < len(texto) and texto[pos] in "-=#":
                enlace = texto[pos]
                pos += 1
            lado.append(Rama(enlace, leer_cadena(hasta=(cierre, ","))))
            saltar_espacios()
            if pos < len(texto) and texto[pos] == ",":
                pos += 1
                saltar_espacios()
                continue
            if pos < len(texto) and texto[pos] == cierre:
                pos += 1
            saltar_espacios()

    def leer_grupo() -> Nodo:
        nonlocal pos
        saltar_espacios()
        inicio = pos
        while pos < len(texto) and (
            texto[pos].isalnum() or texto[pos] in "'⁻⁺"
        ):
            pos += 1
        if inicio == pos:
            raise ValueError(f"Se esperaba un grupo en la posición {pos} de {texto!r}")
        nodo = Nodo(texto[inicio:pos])
        leer_ramas(nodo)
        return nodo

    def leer_cadena(hasta: tuple[str, ...] = ()) -> Cadena:
        nonlocal pos
        cadena: Cadena = []
        while True:
            saltar_espacios()
            if pos >= len(texto) or (hasta and texto[pos] in hasta):
                break
            if texto[pos] in ")]}":
                break
            nodo = leer_grupo()
            saltar_espacios()
            enlace = ""
            if pos < len(texto) and texto[pos] in "-=#":
                enlace = texto[pos]
                pos += 1
            cadena.append((nodo, enlace))
            if not enlace:
                break
        return cadena

    return leer_cadena()


# ---------------------------------------------------------------------------
#  Variantes de la fórmula (para los conmutadores de la aplicación)
# ---------------------------------------------------------------------------

def _expandir_cooh(cadena: Cadena) -> Cadena:
    """Cambia COOH por su forma desarrollada C(=O)-OH.

    El grupo «COOH» se parte en dos: el carbono con su doble enlace al
    oxígeno (dibujado debajo, como en los apuntes) y el OH. El «HOOC» del
    principio (HOOC-CH2-...) se convierte en HO-C(=O)-CH2-.
    """
    salida: Cadena = []
    for nodo, enlace in cadena:
        nodo.arriba = [Rama(r.enlace, _expandir_cooh(r.cadena)) for r in nodo.arriba]
        nodo.abajo = [Rama(r.enlace, _expandir_cooh(r.cadena)) for r in nodo.abajo]
        if not nodo.arriba and not nodo.abajo and nodo.grupo == "COOH":
            carbono = Nodo("C", abajo=[Rama("=", [(Nodo("O"), "")])])
            salida.append((carbono, "-"))
            salida.append((Nodo("OH"), enlace))
        elif not nodo.arriba and not nodo.abajo and nodo.grupo == "HOOC":
            carbono = Nodo("C", abajo=[Rama("=", [(Nodo("O"), "")])])
            salida.append((Nodo("HO"), "-"))
            salida.append((carbono, enlace))
        else:
            salida.append((nodo, enlace))
    return salida


def _es_ch2_simple(nodo: Nodo) -> bool:
    """Un CH2 de la cadena, sin ramas."""
    return nodo.grupo == "CH2" and not nodo.arriba and not nodo.abajo


def _agrupar_ch2(cadena: Cadena) -> Cadena:
    """Junta los CH2 seguidos en un grupo (CH2)n, como en los apuntes.

    CH3-CH2-CH2-CH2-COOH -> CH3-(CH2)3-COOH. Solo se agrupan los CH2
    unidos por enlaces simples (CH2=CH-CH2-CH2-CH3 -> CH2=CH-(CH2)2-CH3).
    """
    salida: Cadena = []
    total = len(cadena)
    i = 0
    while i < total:
        nodo, enlace = cadena[i]
        nodo.arriba = [Rama(r.enlace, _agrupar_ch2(r.cadena)) for r in nodo.arriba]
        nodo.abajo = [Rama(r.enlace, _agrupar_ch2(r.cadena)) for r in nodo.abajo]
        if _es_ch2_simple(nodo):
            j = i
            while (
                j + 1 < total
                and cadena[j][1] == "-"
                and _es_ch2_simple(cadena[j + 1][0])
            ):
                j += 1
            if j > i:
                salida.append((Nodo(f"(CH2){j - i + 1}"), cadena[j][1]))
                i = j + 1
                continue
        salida.append((nodo, enlace))
        i += 1
    return salida


# ---------------------------------------------------------------------------
#  El conmutador «CH3»: fragmentos iguales del mismo átomo, dentro de un grupo
# ---------------------------------------------------------------------------

_GRUPOS_CARBONO = ("C", "CH", "CH2", "CH3")  # nodos de carbono de la cadena


def _es_cadena_simple(nodo: Nodo) -> bool:
    """Un nodo de carbono sin ramas (CH3, CH2, CH)."""
    return (
        nodo.grupo in ("CH3", "CH2", "CH")
        and not nodo.arriba
        and not nodo.abajo
    )


def _tramo_izquierdo(cadena: Cadena, i: int) -> list[int]:
    """Índices del tramo de carbonos simple a la izquierda de `i`, del más
    cercano a `i` hacia fuera."""
    indices: list[int] = []
    j = i - 1
    while j >= 0 and cadena[j][1] == "-" and _es_cadena_simple(cadena[j][0]):
        indices.append(j)
        j -= 1
    return indices


def _tramo_derecho(cadena: Cadena, i: int) -> list[int]:
    """Índices del tramo de carbonos simple a la derecha de `i`, del más
    cercano a `i` hacia fuera."""
    indices: list[int] = []
    j = i + 1
    while (
        j < len(cadena)
        and cadena[j - 1][1] == "-"
        and _es_cadena_simple(cadena[j][0])
    ):
        indices.append(j)
        j += 1
    return indices


def _forma_fragmento(nodos: list[Nodo]) -> str | None:
    """Forma del fragmento escrito desde su extremo libre: «CH3» o «CH3CH2».

    `nodos` va del enlace con el átomo hacia el extremo libre. Devuelve None
    si el fragmento no acaba en CH3 o tiene más de dos carbonos (los apuntes
    solo juntan metilos y etilos).
    """
    if not nodos or len(nodos) > 2 or nodos[-1].grupo != "CH3":
        return None
    if any(n.grupo not in ("CH3", "CH2") for n in nodos):
        return None
    return "".join(n.grupo for n in reversed(nodos))


def _rama_suelta(nodo: Nodo) -> tuple[str, int] | None:
    """La única rama sencilla (un solo grupo, sin ramas) del nodo, si hay
    exactamente una. Son las que se pueden pasar a la cadena (el -NH2 de la
    tert-butilamina, el -CH3 de la dietilmetilamina)."""
    sueltas = [
        (lado, k)
        for lado in ("arriba", "abajo")
        for k, rama in enumerate(getattr(nodo, lado))
        if rama.enlace == "-"
        and len(rama.cadena) == 1
        and not rama.cadena[0][0].arriba
        and not rama.cadena[0][0].abajo
    ]
    return sueltas[0] if len(sueltas) == 1 else None


def _agrupar_ch3(cadena: Cadena) -> Cadena:
    """Junta en un grupo los fragmentos iguales que cuelgan del mismo átomo.

    Como en los apuntes: CH3-CH[CH3]-O-... -> (CH3)2CH-O-...,
    CH3-C{NH2}[CH3]-CH3 -> (CH3)3C-NH2 y
    CH3-CH2-N[CH3]-CH2-CH3 -> (CH3CH2)2N-CH3. Solo se agrupan fragmentos
    que acaban en CH3 y no pasan de dos carbonos: metilos y etilos.

    En un carbono hace falta que algún fragmento agrupado cuelgue de una
    rama: si no, el resultado sería una cadena más larga disfrazada
    (CH3-CH2-CH3 -> (CH3)2CH2, que no se escribe así).
    """
    nodos = list(cadena)
    i = 0
    while i < len(nodos):
        nodo = nodos[i][0]
        nodo.arriba = [Rama(r.enlace, _agrupar_ch3(r.cadena)) for r in nodo.arriba]
        nodo.abajo = [Rama(r.enlace, _agrupar_ch3(r.cadena)) for r in nodo.abajo]

        # Fragmentos que cuelgan del nodo: los tramos de la cadena a cada
        # lado y las ramas. Se guardan por su forma (CH3, CH3CH2...).
        fragmentos: dict[str, list[tuple]] = {}
        izquierda = _tramo_izquierdo(nodos, i)
        if izquierda:
            forma = _forma_fragmento([nodos[j][0] for j in izquierda])
            if forma:
                fragmentos.setdefault(forma, []).append(("izquierda",))
        derecha = _tramo_derecho(nodos, i)
        if derecha:
            forma = _forma_fragmento([nodos[j][0] for j in derecha])
            if forma:
                fragmentos.setdefault(forma, []).append(("derecha",))
        for lado in ("arriba", "abajo"):
            for k, rama in enumerate(getattr(nodo, lado)):
                if rama.enlace != "-":
                    continue
                forma = _forma_fragmento([n for n, _ in rama.cadena])
                if forma:
                    fragmentos.setdefault(forma, []).append((lado, k))

        # La forma con más fragmentos (a igualdad, la más larga).
        candidatas = [
            (forma, recs) for forma, recs in fragmentos.items() if len(recs) >= 2
        ]
        if not candidatas:
            i += 1
            continue
        forma, recs = max(candidatas, key=lambda par: (len(par[1]), len(par[0])))
        if nodo.grupo in _GRUPOS_CARBONO and not any(
            r[0] in ("arriba", "abajo") for r in recs
        ):
            i += 1
            continue

        # Quitar las ramas agrupadas (de la última a la primera).
        for lado in ("arriba", "abajo"):
            indices = sorted((r[1] for r in recs if r[0] == lado), reverse=True)
            for k in indices:
                del getattr(nodo, lado)[k]

        # El grupo se coloca donde estaba un fragmento de la cadena; si todos
        # colgaban de ramas, va pegado a la derecha del nodo. Cuando se
        # agrupan los dos tramos (la dietilmetilamina), se quitan los dos y el
        # grupo ocupa el sitio del izquierdo.
        grupo = Nodo(f"({forma}){len(recs)}")
        hay_izquierda = any(r[0] == "izquierda" for r in recs)
        hay_derecha = any(r[0] == "derecha" for r in recs)
        if hay_izquierda and hay_derecha:
            del nodos[i + 1 : i + 1 + len(derecha)]
            nodos[i] = (nodo, "")
            inicio = i - len(izquierda)
            nodos[inicio:i] = [(grupo, "")]
            i = inicio + 1
        elif hay_izquierda:
            inicio = i - len(izquierda)
            nodos[inicio:i] = [(grupo, "")]
            i = inicio + 1
        elif hay_derecha:
            fin = i + len(derecha)
            enlace_final = nodos[fin][1]
            nodos[i] = (nodo, "")
            nodos[i + 1 : fin + 1] = [(grupo, enlace_final)]
        else:
            nodos.insert(i + 1, (grupo, nodos[i][1]))
            nodos[i] = (nodo, "")

        # Si el átomo se queda sin vecinos en la cadena y solo le cuelga una
        # rama sencilla, la rama pasa a la cadena: (CH3)3C-NH2,
        # (CH3CH2)2N-CH3. Así se escribe en los apuntes.
        if (hay_izquierda and i == len(nodos) - 1) or (hay_derecha and i == 0):
            suelta = _rama_suelta(nodo)
            if suelta is not None:
                lado, k = suelta
                rama = getattr(nodo, lado).pop(k)
                vecino = (rama.cadena[0][0], rama.enlace)
                if hay_izquierda:
                    # pasa al final de la cadena
                    nodos[i] = (nodo, rama.enlace)
                    nodos.insert(i + 1, (vecino[0], ""))
                else:
                    nodos.insert(i, vecino)
                    i += 1
        i += 1
    return nodos


# ---------------------------------------------------------------------------
#  Lienzo: dibuja y va midiendo los límites reales
# ---------------------------------------------------------------------------

class Lienzo:
    """Acumula los elementos del dibujo y calcula la caja que ocupan."""

    def __init__(self) -> None:
        self.elementos: list[str] = []
        self.puntos: list[str] = []
        self.x_min = float("inf")
        self.x_max = float("-inf")
        self.y_min = float("inf")
        self.y_max = float("-inf")

    def _caja(self, x0: float, y0: float, x1: float, y1: float) -> None:
        self.x_min = min(self.x_min, x0, x1)
        self.x_max = max(self.x_max, x0, x1)
        self.y_min = min(self.y_min, y0, y1)
        self.y_max = max(self.y_max, y0, y1)

    def texto(self, fuente: Fuente, texto: str, x: float, y: float, tamano: float) -> float:
        """Dibuja el texto y devuelve el ancho avanzado."""
        ancho = fuente.ancho(texto, tamano)
        self.elementos += fuente.contornos(texto, x, y, tamano)
        self._caja(x, y - tamano * ASCENSO, x + ancho, y + tamano * DESCENSO)
        return ancho

    def linea(self, x1: float, y1: float, x2: float, y2: float) -> None:
        self.elementos.append(
            f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}"/>'
        )
        self._caja(min(x1, x2), min(y1, y2), max(x1, x2), max(y1, y2))

    def puntos_n(self, x: float, y: float, tamano: float) -> None:
        """Los dos puntitos del nitrógeno (el par de electrones)."""
        r = tamano * RADIO_PUNTO
        for signo in (-1, 1):
            cx = x + signo * tamano * SEPARACION_PUNTOS
            self.puntos.append(
                f'<circle cx="{cx:.1f}" cy="{y:.1f}" r="{r:.1f}" style="stroke:none"/>'
            )
            self._caja(cx - r, y - r, cx + r, y + r)

    def svg(self, grosor: float, margen: float) -> str:
        ancho = self.x_max - self.x_min + 2 * margen
        alto = self.y_max - self.y_min + 2 * margen
        dx = margen - self.x_min
        dy = margen - self.y_min
        cuerpo = "\n".join(self.elementos)
        puntos = ""
        if self.puntos:
            # van en un grupo aparte: la aplicación los enseña o los esconde
            puntos = (
                f'<g transform="translate({dx:.1f},{dy:.1f})" class="puntos-n" '
                f'fill="#000000">' + "".join(self.puntos) + "</g>"
            )
        return (
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{ancho:.0f}" '
            f'height="{alto:.0f}" viewBox="0 0 {ancho:.0f} {alto:.0f}">'
            f'<rect width="100%" height="100%" fill="#ffffff"/>'
            f'<g transform="translate({dx:.1f},{dy:.1f})" stroke="#000000" '
            f'stroke-width="{grosor:.1f}" stroke-linecap="round" fill="#000000">'
            f"{cuerpo}</g>{puntos}</svg>"
        )


# ---------------------------------------------------------------------------
#  Medidas
# ---------------------------------------------------------------------------

@dataclass
class Medidas:
    ancho: float
    arriba: float
    abajo: float


@dataclass
class MedidasNodo:
    izq: float    # avance a la izquierda del centro (solo el texto)
    der: float    # avance a la derecha del centro (solo el texto)
    arriba: float
    abajo: float
    ramas_arriba_izq: float = 0.0
    ramas_arriba_der: float = 0.0
    ramas_abajo_izq: float = 0.0
    ramas_abajo_der: float = 0.0


def _ancho_grupo(fuente: Fuente, nodo: Nodo, tamano: float) -> float:
    total = 0.0
    for trozo in nodo.trozos():
        if trozo.texto == "⁻":
            total += tamano * ANCHO_SUPER
            continue
        total += fuente.ancho(trozo.texto, tamano * (FS_SUB if trozo.subindice else 1.0))
    return total


def _primera_parte(grupo: str) -> str:
    """Primer átomo del grupo: «Cl»/«Br» van juntos, el resto es una letra."""
    return grupo[:2] if grupo[:2] in ("Cl", "Br") else grupo[:1]


def _ancla_grupo(fuente: Fuente, nodo: Nodo, tamano: float) -> float:
    """Centro del primer átomo respecto al centro del grupo.

    Los enlaces verticales salen del primer átomo (la C del grupo), no del
    centro del texto: así la línea empieza y acaba en la C, como en los apuntes.
    """
    ancho = _ancho_grupo(fuente, nodo, tamano)
    primera = fuente.ancho(_primera_parte(nodo.grupo), tamano)
    return primera / 2 - ancho / 2


def _medir_nodo(fuente: Fuente, nodo: Nodo, tamano: float) -> MedidasNodo:
    """Mide un nodo. Las ramas no ensanchan el hueco del grupo: van por encima
    y por debajo, así que el texto de al lado puede pasar por debajo de ellas.
    Solo se guarda cuánto sobresalen por cada lado, para separar dos ramas
    vecinas que sí chocarían."""
    ancho_texto = _ancho_grupo(fuente, nodo, tamano)
    m = MedidasNodo(ancho_texto / 2, ancho_texto / 2, tamano * ASCENSO, tamano * DESCENSO)
    ancla = _ancla_grupo(fuente, nodo, tamano)

    for lado, ramas in ((True, nodo.arriba), (False, nodo.abajo)):
        for rama in ramas:
            if not rama.cadena:
                continue
            # la primera letra de la rama queda centrada bajo el primer átomo
            # del grupo padre; se mide cuánto sobresale la rama por cada lado
            primera = fuente.ancho(_primera_parte(rama.cadena[0][0].grupo), tamano) / 2
            rm = medir(fuente, rama.cadena, tamano)
            izq = primera - ancla
            der = max(0.0, rm.ancho - izq)
            if lado:
                m.ramas_arriba_izq = max(m.ramas_arriba_izq, max(0.0, izq))
                m.ramas_arriba_der = max(m.ramas_arriba_der, der)
                m.arriba = max(m.arriba, tamano * RAMAS_ARRIBA + rm.arriba)
            else:
                m.ramas_abajo_izq = max(m.ramas_abajo_izq, max(0.0, izq))
                m.ramas_abajo_der = max(m.ramas_abajo_der, der)
                m.abajo = max(m.abajo, tamano * RAMAS_ABAJO + rm.abajo)

    return m


def _disponer(
    fuente: Fuente, cadena: Cadena, tamano: float
) -> tuple[list[MedidasNodo], list[float], float]:
    """Coloca los nodos de la cadena.

    Devuelve las medidas de cada nodo, la posición del centro de cada uno
    (relativa al inicio) y el ancho total. Si dos ramas seguidas van del mismo
    lado y se solaparían, se añade el hueco necesario entre los dos grupos.
    """
    medidas = [_medir_nodo(fuente, nodo, tamano) for nodo, _ in cadena]
    extras: list[float] = []
    for i in range(len(cadena) - 1):
        m1, m2 = medidas[i], medidas[i + 1]
        enlace = cadena[i][1]
        hueco = m1.der + (tamano * ENLACE if enlace else 0.0) + m2.izq
        # Solo hay choque si los dos grupos tienen rama del mismo lado:
        # el texto de al lado va a la altura de la cadena, no de la rama.
        # Se deja además un hueco mínimo entre las dos ramas.
        necesario = 0.0
        hueco_ramas = tamano * HUECO_RAMAS
        if m1.ramas_arriba_der and m2.ramas_arriba_izq:
            necesario = max(necesario, m1.ramas_arriba_der + m2.ramas_arriba_izq + hueco_ramas)
        if m1.ramas_abajo_der and m2.ramas_abajo_izq:
            necesario = max(necesario, m1.ramas_abajo_der + m2.ramas_abajo_izq + hueco_ramas)
        extras.append(max(0.0, necesario - hueco))

    centros: list[float] = []
    cursor = 0.0
    for i, (nodo, enlace) in enumerate(cadena):
        m = medidas[i]
        centros.append(cursor + m.izq)
        cursor += m.izq + m.der
        if enlace and i < len(cadena) - 1:
            cursor += tamano * ENLACE + extras[i]
    return medidas, centros, cursor


def centros_grupos(fuente: Fuente, cadena: Cadena, tamano: float = FS) -> list[float]:
    """Centro (x) de cada grupo de la cadena, relativo al inicio."""
    _, centros, _ = _disponer(fuente, cadena, tamano)
    return centros


def medir(fuente: Fuente, cadena: Cadena, tamano: float = FS) -> Medidas:
    if not cadena:
        return Medidas(0.0, tamano * ASCENSO, tamano * DESCENSO)
    medidas, _, ancho = _disponer(fuente, cadena, tamano)
    arriba = max([tamano * ASCENSO] + [m.arriba for m in medidas])
    abajo = max([tamano * DESCENSO] + [m.abajo for m in medidas])
    return Medidas(ancho, arriba, abajo)


# ---------------------------------------------------------------------------
#  Dibujo
# ---------------------------------------------------------------------------

def _enlace_horizontal(
    lienzo: Lienzo, enlace: str, x1: float, x2: float, y: float, tamano: float
) -> None:
    if enlace == "-":
        lienzo.linea(x1, y, x2, y)
    elif enlace == "=":
        # el doble va centrado en la altura de las letras (como el triple),
        # con las dos rayitas a la misma distancia
        centro = y + tamano * (Y_ENLACE - Y_TRIPLE)
        for d in (-1, 1):
            lienzo.linea(
                x1, centro + d * tamano * SEPARACION_DOBLE,
                x2, centro + d * tamano * SEPARACION_DOBLE,
            )
    else:
        # el triple va centrado en la altura de las letras, no en el eje
        centro = y + tamano * (Y_ENLACE - Y_TRIPLE)
        for d in (-1, 0, 1):
            lienzo.linea(
                x1, centro + d * tamano * SEPARACION_TRIPLE,
                x2, centro + d * tamano * SEPARACION_TRIPLE,
            )


def _enlace_vertical(
    lienzo: Lienzo, enlace: str, x: float, y1: float, y2: float, tamano: float
) -> None:
    if enlace == "-":
        lienzo.linea(x, y1, x, y2)
    elif enlace == "=":
        for d in (-1, 1):
            lienzo.linea(x + d * tamano * SEPARACION_DOBLE, y1, x + d * tamano * SEPARACION_DOBLE, y2)
    else:
        centro = x + tamano * (Y_ENLACE - Y_TRIPLE)
        for d in (-1, 0, 1):
            lienzo.linea(
                centro + d * tamano * SEPARACION_TRIPLE, y1,
                centro + d * tamano * SEPARACION_TRIPLE, y2,
            )


def _texto_nodo(
    lienzo: Lienzo, fuente: Fuente, nodo: Nodo, centro: float, base: float, tamano: float
) -> None:
    ancho = _ancho_grupo(fuente, nodo, tamano)
    cursor = centro - ancho / 2
    x_nitrogeno = None
    for trozo in nodo.trozos():
        if trozo.texto == "⁻":
            # la rayita del superíndice (O-), que la fuente no tiene
            fin = cursor + tamano * ANCHO_SUPER
            y = base - tamano * SUBE_SUPER
            lienzo.linea(cursor + tamano * 0.04, y, fin - tamano * 0.04, y)
            cursor = fin
            continue
        t = tamano * (FS_SUB if trozo.subindice else 1.0)
        y = base + (tamano * SUB_BAJADA if trozo.subindice else 0.0)
        if (
            x_nitrogeno is None
            and not trozo.subindice
            and "N" in trozo.texto
            and not nodo.grupo.startswith("NO")  # el nitro no lleva puntitos
        ):
            antes = trozo.texto[: trozo.texto.index("N")]
            x_nitrogeno = cursor + fuente.ancho(antes, t) + fuente.ancho("N", t) / 2
        cursor += lienzo.texto(fuente, trozo.texto, cursor, y, t)
    if x_nitrogeno is not None:
        lienzo.puntos_n(x_nitrogeno, base - tamano * (ALTO_GLIFO + 0.12), tamano)


def _dibujar_nodo(
    lienzo: Lienzo, fuente: Fuente, nodo: Nodo, centro: float, base: float, tamano: float
) -> None:
    """Dibuja un grupo con sus ramas. `centro` es el centro del grupo."""
    _texto_nodo(lienzo, fuente, nodo, centro, base, tamano)
    ancla = centro + _ancla_grupo(fuente, nodo, tamano)

    for rama in nodo.arriba:
        if not rama.cadena:
            continue
        base_rama = base - tamano * RAMAS_ARRIBA
        _enlace_vertical(
            lienzo, rama.enlace, ancla,
            base - tamano * INICIO_ARRIBA, base_rama + tamano * FIN_RAMA, tamano,
        )
        primera = fuente.ancho(_primera_parte(rama.cadena[0][0].grupo), tamano) / 2
        dibujar(lienzo, fuente, rama.cadena, ancla - primera, base_rama, tamano)

    for rama in nodo.abajo:
        if not rama.cadena:
            continue
        base_rama = base + tamano * RAMAS_ABAJO
        _enlace_vertical(
            lienzo, rama.enlace, ancla,
            base + tamano * INICIO_ABAJO,
            base_rama - tamano * (ALTO_GLIFO + FIN_RAMA), tamano,
        )
        primera = fuente.ancho(_primera_parte(rama.cadena[0][0].grupo), tamano) / 2
        dibujar(lienzo, fuente, rama.cadena, ancla - primera, base_rama, tamano)


def dibujar(
    lienzo: Lienzo, fuente: Fuente, cadena: Cadena, x: float, base: float, tamano: float
) -> None:
    """Dibuja una cadena horizontal: `x` es el borde izquierdo del primer grupo."""
    if not cadena:
        return
    medidas, centros, _ = _disponer(fuente, cadena, tamano)
    for i, (nodo, enlace) in enumerate(cadena):
        centro = x + centros[i]
        _dibujar_nodo(lienzo, fuente, nodo, centro, base, tamano)
        if enlace and i < len(cadena) - 1:
            m1, m2 = medidas[i], medidas[i + 1]
            ini = x + centros[i] + m1.der
            fin = x + centros[i + 1] - m2.izq
            margen = (fin - ini) * 0.12
            _enlace_horizontal(
                lienzo, enlace, ini + margen, fin - margen,
                base - tamano * Y_ENLACE, tamano,
            )


def svg_condensada(
    formula: str,
    tamano: float = FS,
    *,
    cooh_separado: bool = False,
    ch2_agrupados: bool = False,
    ch3_agrupados: bool = False,
) -> str:
    """Devuelve el SVG de la fórmula semidesarrollada.

    `cooh_separado` dibuja los COOH como C(=O)-OH, `ch2_agrupados` junta los
    CH2 seguidos en (CH2)n y `ch3_agrupados` junta en un grupo los CH3 (o los
    CH3CH2) que cuelgan del mismo átomo: son los conmutadores de la
    aplicación.
    """
    fuente = Fuente()
    formula = formula.strip()
    abierto_izq = formula.startswith("-")
    abierto_der = formula.endswith("-")
    formula = formula.strip("-").strip()
    cadena = analizar(formula)
    if cooh_separado:
        cadena = _expandir_cooh(cadena)
    if ch2_agrupados:
        cadena = _agrupar_ch2(cadena)
    if ch3_agrupados:
        cadena = _agrupar_ch3(cadena)

    lienzo = Lienzo()
    hueco = tamano * ENLACE * 0.85
    x = 0.0
    if abierto_izq:
        _enlace_horizontal(lienzo, "-", 0.0, hueco, -tamano * Y_ENLACE, tamano)
        x += hueco
    dibujar(lienzo, fuente, cadena, x, 0.0, tamano)
    if abierto_der:
        fin = lienzo.x_max
        _enlace_horizontal(lienzo, "-", fin + hueco * 0.15, fin + hueco * 1.0,
                           -tamano * Y_ENLACE, tamano)

    return lienzo.svg(grosor=tamano * GROSOR, margen=tamano * MARGEN)


if __name__ == "__main__":
    salida = Path("/tmp/opencode/salida_condensada")
    salida.mkdir(parents=True, exist_ok=True)
    pruebas = {
        "ejemplo1": "CH3-CH{F}-CH2-F",
        "butanona": "CH3-C{=O}-CH2-CH3",
        "dimetilpentanona": "CH3-C{=O}-CH2-C{CH3}[CH3]-CH3",
        "metilbuteno": "CH3-C{CH3}=CH-CH3",
        "hepteno": "CH3-CH=C{CH3}-CH2-CH2-CH{CH3}-CH3",
        "heptino": "CH3-CH{CH3}-C#C-CH2-CH{Br}-CH3",
        "pentenol": "CH3-CH2-C{OH}=CH-CH3",
        "vinilo": "-CH=CH2",
        "cetona": "R-C{=O}-R'",
        "butanol": "CH3-C{CH3}[CH3]-CH{OH}-CH2-Br",
        "propilhepteno": "CH2=CH-CH{CH2-CH2-CH3}-CH2-CH2-CH2-CH3",
        "butadieno": "CH2=CH-CH=CH2",
        "etino": "CH#CH",
        "acidooxobutanoico": "CH3-C{=O}-CH2-COOH",
        "dimetilhepteno": "CH3-CH=C{CH3}-CH2-CH2-CH{CH3}-CH3",
    }
    for nombre, formula in pruebas.items():
        svg = svg_condensada(formula)
        (salida / f"{nombre}.svg").write_text(svg, encoding="utf-8")
        print(f"{nombre}: {len(svg)} bytes")
