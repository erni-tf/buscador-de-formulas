# Buscador de Fórmulas

La química orgánica de los apuntes de clase, buscable y dibujada como en el
papel. Este glosario fija el vocabulario con el que se habla del catálogo, de
las vistas y de la postura de cada molécula.

## Lenguaje

### El catálogo y su origen

**Apuntes**:
Las fotos de clase de las que salen los compuestos, los grupos y la forma en
que se dibuja cada uno. Son la fuente de verdad del proyecto.
_Avoid_: notas, temario, material

**Catálogo**:
La base de datos editable de entradas (`datos/compuestos.toml`).
_Avoid_: base de datos, dataset, inventario

**Entrada**:
Cada bloque del catálogo: un compuesto o una ficha.
_Avoid_: registro, elemento, item

**Compuesto**:
Entrada con una estructura química propia, buscable por nombre, alias o
fórmula molecular.
_Avoid_: molécula (cuando se habla de la entrada)

**Ficha**:
Entrada que enseña un grupo funcional o un concepto genérico, con comodines
en lugar de una molécula concreta.
_Avoid_: tarjeta

**Familia**:
Agrupación del catálogo por tipo de compuesto (alcanos, alcoholes…). Una
entrada puede estar en varias.
_Avoid_: categoría, grupo (para la agrupación del filtro)

### Las vistas

**Vista**:
Cada uno de los dibujos de una entrada: la estructura y la fórmula
semidesarrollada.
_Avoid_: pestaña, modo

**Estructura**:
La vista de líneas, con los heteroátomos escritos y los carbonos implícitos.
_Avoid_: esqueleto (cuando se habla de la vista), dibujo de líneas

**Fórmula semidesarrollada**:
La vista de texto con los grupos escritos (`CH3-C(=O)-CH2-CH3`), dibujada como
en los apuntes.
_Avoid_: fórmula condensada, condensada

**Conmutador**:
Control de pantalla que cambia la forma de una vista: el anillo aromático con
rayas o con círculo, el par de electrones del nitrógeno, los COOH juntos o
desarrollados y los CH₂ seguidos o los CH₃ (y etilos) repetidos sueltos o
agrupados.
_Avoid_: opción, ajuste, toggle

### La postura

**Postura**:
La colocación de los átomos de una molécula antes de dibujarla, para que el
dibujo quede como en los apuntes.
_Avoid_: orientación, disposición, layout

**Regla de postura**:
Una forma nombrada de colocar la molécula, con la condición que dice a qué
moléculas se aplica. Cada regla reconoce sus moléculas y las coloca; la
primera que las reclama es la que actúa.
_Avoid_: transformación, heurística, caso especial

**Postura por defecto**:
La colocación que deja RDKit cuando ninguna regla de postura reclama la
molécula.
_Avoid_: postura automática, sin regla

**Quiebro**:
El zigzag que se da a las cadenas dibujadas rectas (los alquinos), para que
cada enlace cambie de dirección.
_Avoid_: doblado, ángulo

**Espejo**:
El reflejo horizontal del dibujo, para que un sustituyente escrito (el
vinilo) caiga a la derecha de su unión.
_Avoid_: volteo, inversión
