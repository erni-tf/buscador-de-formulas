/* ==========================================================================
   Buscador de Fórmulas — lógica de la aplicación
   ========================================================================== */
(function () {
  "use strict";

  var datos = window.DATOS || { entradas: [], familias: [] };
  var entradas = datos.entradas;

  var elQ = document.getElementById("q");
  var elLimpiar = document.getElementById("limpiar");
  var elSug = document.getElementById("sugerencias");
  var elFiltros = document.getElementById("filtros");
  var elDetalle = document.getElementById("detalle");
  var elLista = document.getElementById("lista");
  var elVacio = document.getElementById("vacio");
  var elNoEncontrado = document.getElementById("no-encontrado");
  var elTituloLista = document.getElementById("titulo-lista");

  var familiaActiva = null;
  var sugerencias = [];
  var sugerenciaActiva = -1;

  /* Familias de una entrada (una entrada puede salir en varias páginas). */
  function familiasDe(e) {
    return e.familias && e.familias.length ? e.familias : [e.familia];
  }

  /* Anillos aromáticos: con dobles enlaces alternos («3 rayas») o con
     círculo. La elección se recuerda entre visitas. */
  var circulos = false;
  try {
    circulos = window.localStorage.getItem("aromaticos") === "circulo";
  } catch (err) { circulos = false; }

  function guardarCirculos() {
    try {
      window.localStorage.setItem("aromaticos", circulos ? "circulo" : "rayas");
    } catch (err) { /* si no se puede guardar, no pasa nada */ }
  }

  /* Los dos puntitos del nitrógeno (el par de electrones). Por defecto se
     enseñan, como en los apuntes. */
  var puntosN = true;
  try {
    puntosN = window.localStorage.getItem("puntosN") !== "sin";
  } catch (err) { puntosN = true; }

  function guardarPuntosN() {
    try {
      window.localStorage.setItem("puntosN", puntosN ? "con" : "sin");
    } catch (err) { /* si no se puede guardar, no pasa nada */ }
  }

  /* Los COOH de las fórmulas semidesarrolladas: juntos (COOH) o separados
     (C(=O)-OH). Por defecto, juntos, como en los apuntes. */
  var coohSeparado = false;
  try {
    coohSeparado = window.localStorage.getItem("cooh") === "separado";
  } catch (err) { coohSeparado = false; }

  function guardarCooh() {
    try {
      window.localStorage.setItem("cooh", coohSeparado ? "separado" : "junto");
    } catch (err) { /* si no se puede guardar, no pasa nada */ }
  }

  /* Los CH2 seguidos: sueltos (CH2-CH2-CH2) o juntos en un paréntesis
     ((CH2)3), como en los apuntes del ácido pentanodioico. */
  var ch2Agrupados = false;
  try {
    ch2Agrupados = window.localStorage.getItem("ch2") === "agrupados";
  } catch (err) { ch2Agrupados = false; }

  function guardarCh2() {
    try {
      window.localStorage.setItem("ch2", ch2Agrupados ? "agrupados" : "sueltos");
    } catch (err) { /* si no se puede guardar, no pasa nada */ }
  }

  /* ---------------------------------------------------------------- texto */

  function sinAcentos(texto) {
    return String(texto).toLowerCase()
      .replace(/[áàäâã]/g, "a").replace(/[éèëê]/g, "e")
      .replace(/[íìïî]/g, "i").replace(/[óòöô]/g, "o")
      .replace(/[úùüû]/g, "u").replace(/ñ/g, "n");
  }

  function clave(texto) {
    return sinAcentos(texto).replace(/[^a-z0-9]/g, "");
  }

  /* Fórmulas bonitas: los dígitos como subíndices (C6H10 -> C₆H₁₀). */
  var SUBINDICES = {
    "0": "₀", "1": "₁", "2": "₂", "3": "₃", "4": "₄",
    "5": "₅", "6": "₆", "7": "₇", "8": "₈", "9": "₉"
  };

  function conSubindices(texto) {
    return String(texto == null ? "" : texto).replace(/\d/g, function (d) {
      return SUBINDICES[d];
    });
  }

  function distancia(a, b) {
    if (Math.abs(a.length - b.length) > 2) return 99;
    var fila = [];
    var i, j;
    for (j = 0; j <= b.length; j++) fila[j] = j;
    for (i = 1; i <= a.length; i++) {
      var anterior = fila[0];
      fila[0] = i;
      for (j = 1; j <= b.length; j++) {
        var temp = fila[j];
        fila[j] = Math.min(
          fila[j] + 1,
          fila[j - 1] + 1,
          anterior + (a.charAt(i - 1) === b.charAt(j - 1) ? 0 : 1)
        );
        anterior = temp;
      }
    }
    return fila[b.length];
  }

  /* --------------------------------------------------------- preparativos */

  entradas.forEach(function (e) {
    var nombres = e.nombres && e.nombres.length ? e.nombres : [e.nombre];
    e._claves = nombres.map(clave);
    e._palabras = nombres.map(function (n) {
      return sinAcentos(n).split(/[\s\-]+/).map(clave).filter(function (p) { return p; });
    });
    e._formulaClave = clave(e.formula);
  });

  function imagen(entrada, tipo) {
    var nodo = document.getElementById("img-" + tipo + "-" + entrada.id);
    return nodo ? nodo.innerHTML : "";
  }

  /* ------------------------------------------------------------- búsqueda */

  function buscar(consulta) {
    var cq = clave(consulta);
    if (!cq) return [];
    var res = [];
    for (var i = 0; i < entradas.length; i++) {
      var e = entradas[i];
      var mejor = 99;
      for (var j = 0; j < e._claves.length; j++) {
        var k = e._claves[j];
        if (k === cq) { mejor = 0; break; }
        if (k.indexOf(cq) === 0) { mejor = Math.min(mejor, 1); continue; }
        var palabras = e._palabras[j];
        var porPalabra = false;
        for (var w = 0; w < palabras.length; w++) {
          if (palabras[w].indexOf(cq) === 0) { porPalabra = true; break; }
        }
        if (porPalabra) { mejor = Math.min(mejor, 2); continue; }
        if (k.indexOf(cq) > 0) { mejor = Math.min(mejor, 3); continue; }
        if (cq.length >= 4 && distancia(k, cq) <= 2) mejor = Math.min(mejor, 4);
      }
      if (e._formulaClave && e._formulaClave === cq) mejor = 0;
      if (mejor < 99) res.push({ e: e, p: mejor });
    }
    res.sort(function (a, b) {
      if (a.p !== b.p) return a.p - b.p;
      var la = a.e.nombre.length, lb = b.e.nombre.length;
      if (la !== lb) return la - lb;
      return a.e.nombre.localeCompare(b.e.nombre, "es");
    });
    return res;
  }

  /* -------------------------------------------------------------- listado */

  function visibles() {
    if (!familiaActiva) return entradas;
    return entradas.filter(function (e) { return familiasDe(e).indexOf(familiaActiva) >= 0; });
  }

  function pintarFiltros() {
    elFiltros.innerHTML = "";
    datos.familias.forEach(function (f) {
      var b = document.createElement("button");
      b.type = "button";
      b.className = "chip" + (familiaActiva === f ? " activo" : "");
      b.textContent = f;
      b.onclick = function () {
        familiaActiva = (familiaActiva === f) ? null : f;
        pintarFiltros();
        pintarLista();
      };
      elFiltros.appendChild(b);
    });
  }

  function pintarLista() {
    var lista = visibles();
    elTituloLista.textContent = familiaActiva
      ? "Familia: " + familiaActiva + " (" + lista.length + ")"
      : "Todos los compuestos (" + lista.length + ")";
    elLista.innerHTML = "";
    lista.forEach(function (e) {
      var b = document.createElement("button");
      b.type = "button";
      b.className = "item";
      var n = document.createElement("span");
      n.textContent = e.nombre;
      b.appendChild(n);
      var f = document.createElement("span");
      f.className = "fam";
      f.textContent = familiasDe(e).join(" · ") + (e.formula ? " · " + conSubindices(e.formula) : "");
      b.appendChild(f);
      b.onclick = function () { mostrar(e); };
      elLista.appendChild(b);
    });
  }

  /* -------------------------------------------------------------- detalle */

  function mostrar(e) {
    elNoEncontrado.hidden = true;
    elVacio.hidden = true;
    elDetalle.hidden = false;
    elDetalle.innerHTML = "";
    elDetalle.className = "tarjeta";

    var cab = document.createElement("div");
    cab.className = "detalle-cabecera";
    var izq = document.createElement("div");
    var h = document.createElement("h2");
    h.className = "detalle-nombre";
    h.textContent = e.nombre;
    izq.appendChild(h);
    if (e.alias && e.alias.length) {
      var al = document.createElement("p");
      al.className = "detalle-alias";
      al.textContent = "También: " + e.alias.join(" · ");
      izq.appendChild(al);
    }
    var meta = document.createElement("div");
    meta.className = "detalle-meta";
    var et = document.createElement("span");
    et.className = "etiqueta";
    et.textContent = familiasDe(e).join(" · ");
    meta.appendChild(et);
    if (e.formula) {
      var fm = document.createElement("span");
      fm.className = "formula-molecular";
      fm.textContent = conSubindices(e.formula);
      meta.appendChild(fm);
    }
    if (e.tipo === "ficha") {
      var fi = document.createElement("span");
      fi.textContent = "ficha de " + familiasDe(e).join(" · ");
      meta.appendChild(fi);
    }
    izq.appendChild(meta);
    cab.appendChild(izq);
    elDetalle.appendChild(cab);

    var vistas = document.createElement("div");
    vistas.className = "vistas";

    function añadirVista(titulo, contenido, tipo) {
      var v = document.createElement("div");
      v.className = "vista";
      var t = document.createElement("div");
      t.className = "titulo";
      t.textContent = titulo;
      v.appendChild(t);
      var zona = document.createElement("div");
      zona.className = "zona";
      zona.innerHTML = contenido;
      v.appendChild(zona);
      vistas.appendChild(v);
      aplicarPuntosN(zona);
      var svg = v.querySelector("svg");
      if (svg) v.appendChild(accionesImagen(svg, nombreArchivo(e.id, tipo)));
      return v;
    }

    /* Los puntitos de la N se enseñan o se esconden (también al copiar) */
    function aplicarPuntosN(raiz) {
      var grupos = raiz.querySelectorAll(".puntos-n");
      for (var k = 0; k < grupos.length; k++) {
        grupos[k].style.display = puntosN ? "" : "none";
      }
    }

    /* Los dibujos que hay que enseñar los resuelve el build: aquí sólo se
       lee la casilla de la combinación de conmutadores en la que estamos. */
    function tipoEstructura() {
      return e.vistas.estructura[circulos ? "circulo" : "rayas"];
    }

    function combinacion() {
      return (circulos ? "circulo" : "rayas") + "-" +
        (coohSeparado ? "separado" : "junto") + "-" +
        (ch2Agrupados ? "agrupados" : "sueltos");
    }

    function pintarVistas() {
      vistas.innerHTML = "";
      añadirVista("Estructura", imagen(e, tipoEstructura()), "estructura");
      if (e.vistas.semidesarrollada) {
        añadirVista(
          "Fórmula semidesarrollada",
          imagen(e, e.vistas.semidesarrollada[combinacion()]),
          "semidesarrollada"
        );
      }
    }

    if (e.circulos) {
      var conm = document.createElement("div");
      conm.className = "conmutador";
      var ct = document.createElement("span");
      ct.className = "conmutador-titulo";
      ct.textContent = "Anillo aromático:";
      conm.appendChild(ct);
      [["rayas", "3 rayas"], ["circulo", "círculo"]].forEach(function (op) {
        var b = document.createElement("button");
        b.type = "button";
        b.className = "chip" + ((circulos ? "circulo" : "rayas") === op[0] ? " activo" : "");
        b.textContent = op[1];
        b.onclick = function () {
          circulos = (op[0] === "circulo");
          guardarCirculos();
          mostrar(e);
        };
        conm.appendChild(b);
      });
      elDetalle.appendChild(conm);
    }

    if (e.puntos) {
      var conmN = document.createElement("div");
      conmN.className = "conmutador";
      var ctN = document.createElement("span");
      ctN.className = "conmutador-titulo";
      ctN.textContent = "Nitrógeno:";
      conmN.appendChild(ctN);
      [["puntos", "con 2 puntos"], ["sin", "sin 2 puntos"]].forEach(function (op) {
        var b = document.createElement("button");
        b.type = "button";
        b.className = "chip" + ((puntosN ? "puntos" : "sin") === op[0] ? " activo" : "");
        b.textContent = op[1];
        b.onclick = function () {
          puntosN = (op[0] === "puntos");
          guardarPuntosN();
          mostrar(e);
        };
        conmN.appendChild(b);
      });
      elDetalle.appendChild(conmN);
    }

    if (e.cooh) {
      var conmC = document.createElement("div");
      conmC.className = "conmutador";
      var ctC = document.createElement("span");
      ctC.className = "conmutador-titulo";
      ctC.textContent = "Grupo COOH:";
      conmC.appendChild(ctC);
      [["junto", "COOH"], ["separado", "C(=O)-OH"]].forEach(function (op) {
        var b = document.createElement("button");
        b.type = "button";
        b.className = "chip" + ((coohSeparado ? "separado" : "junto") === op[0] ? " activo" : "");
        b.textContent = op[1];
        b.onclick = function () {
          coohSeparado = (op[0] === "separado");
          guardarCooh();
          mostrar(e);
        };
        conmC.appendChild(b);
      });
      elDetalle.appendChild(conmC);
    }

    if (e.ch2) {
      var conmCh = document.createElement("div");
      conmCh.className = "conmutador";
      var ctCh = document.createElement("span");
      ctCh.className = "conmutador-titulo";
      ctCh.textContent = "CH₂ seguidos:";
      conmCh.appendChild(ctCh);
      [["sueltos", "sueltos"], ["agrupados", "(CH₂)ₙ"]].forEach(function (op) {
        var b = document.createElement("button");
        b.type = "button";
        b.className = "chip" + ((ch2Agrupados ? "agrupados" : "sueltos") === op[0] ? " activo" : "");
        b.textContent = op[1];
        b.onclick = function () {
          ch2Agrupados = (op[0] === "agrupados");
          guardarCh2();
          mostrar(e);
        };
        conmCh.appendChild(b);
      });
      elDetalle.appendChild(conmCh);
    }

    pintarVistas();
    elDetalle.appendChild(vistas);

    if (e.notas) {
      var n = document.createElement("div");
      n.className = "notas";
      n.textContent = e.notas;
      elDetalle.appendChild(n);
    }

    if (e.condensada_texto) {
      var bts = document.createElement("div");
      bts.className = "botones";
      var bc = document.createElement("button");
      bc.type = "button";
      bc.className = "boton";
      bc.textContent = "Copiar la fórmula";
      bc.onclick = function () {
        copiar(e.condensada_texto);
        bc.textContent = "¡Copiada!";
        setTimeout(function () { bc.textContent = "Copiar la fórmula"; }, 1400);
      };
      bts.appendChild(bc);
      elDetalle.appendChild(bts);
    }

    if (window.scrollY > 10) {
      window.scrollTo({ top: 0, behavior: "smooth" });
    }
  }

  function copiar(texto) {
    try {
      if (navigator.clipboard && navigator.clipboard.writeText) {
        navigator.clipboard.writeText(texto);
        return;
      }
    } catch (err) { /* seguimos con el método clásico */ }
    var ta = document.createElement("textarea");
    ta.value = texto;
    ta.style.position = "fixed";
    ta.style.opacity = "0";
    document.body.appendChild(ta);
    ta.select();
    try { document.execCommand("copy"); } catch (err2) { /* nada */ }
    document.body.removeChild(ta);
  }

  /* --------------------------------------------------- copiar la imagen */

  function svgABlob(svg, escala, alTerminar) {
    var ancho = parseFloat(svg.getAttribute("width")) || 0;
    var alto = parseFloat(svg.getAttribute("height")) || 0;
    if (!ancho || !alto) {
      try {
        ancho = svg.viewBox.baseVal.width;
        alto = svg.viewBox.baseVal.height;
      } catch (err) { /* nos quedamos con 0 */ }
    }
    if (!ancho || !alto) { ancho = 440; alto = 320; }

    var datos = new XMLSerializer().serializeToString(svg);
    if (datos.indexOf("xmlns=") === -1) {
      datos = datos.replace("<svg", '<svg xmlns="http://www.w3.org/2000/svg"');
    }
    var url;
    try {
      url = URL.createObjectURL(new Blob([datos], { type: "image/svg+xml;charset=utf-8" }));
    } catch (err) {
      alTerminar(null);
      return;
    }
    var img = new Image();
    img.onload = function () {
      try {
        var c = document.createElement("canvas");
        c.width = Math.round(ancho * escala);
        c.height = Math.round(alto * escala);
        var ctx = c.getContext("2d");
        ctx.fillStyle = "#ffffff";
        ctx.fillRect(0, 0, c.width, c.height);
        ctx.drawImage(img, 0, 0, c.width, c.height);
        URL.revokeObjectURL(url);
        if (c.toBlob) {
          c.toBlob(function (b) { alTerminar(b); }, "image/png");
        } else {
          alTerminar(null);
        }
      } catch (err) {
        alTerminar(null);
      }
    };
    img.onerror = function () {
      URL.revokeObjectURL(url);
      alTerminar(null);
    };
    img.src = url;
  }

  function descargar(blob, nombre) {
    var url = URL.createObjectURL(blob);
    var a = document.createElement("a");
    a.href = url;
    a.download = nombre;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    setTimeout(function () { URL.revokeObjectURL(url); }, 4000);
  }

  function aviso(boton, texto) {
    var original = boton.getAttribute("data-texto") || boton.textContent;
    boton.setAttribute("data-texto", original);
    boton.textContent = texto;
    setTimeout(function () { boton.textContent = original; }, 1600);
  }

  function accionesImagen(svg, nombre) {
    var caja = document.createElement("div");
    caja.className = "acciones-imagen";

    var bc = document.createElement("button");
    bc.type = "button";
    bc.className = "enlace";
    bc.textContent = "Copiar imagen";
    bc.onclick = function () {
      svgABlob(svg, 3, function (blob) {
        if (!blob) { aviso(bc, "No se pudo"); return; }
        var copiado = false;
        try {
          if (navigator.clipboard && window.ClipboardItem) {
            navigator.clipboard.write([new ClipboardItem({ "image/png": blob })])
              .then(function () { aviso(bc, "¡Copiada!"); })
              .catch(function () { descargar(blob, nombre); aviso(bc, "Descargada"); });
            copiado = true;
          }
        } catch (err) { copiado = false; }
        if (!copiado) {
          descargar(blob, nombre);
          aviso(bc, "Descargada");
        }
      });
    };
    caja.appendChild(bc);

    var sep = document.createElement("span");
    sep.className = "sep";
    sep.textContent = "·";
    caja.appendChild(sep);

    var bd = document.createElement("button");
    bd.type = "button";
    bd.className = "enlace";
    bd.textContent = "Descargar PNG";
    bd.onclick = function () {
      svgABlob(svg, 3, function (blob) {
        if (!blob) { aviso(bd, "No se pudo"); return; }
        descargar(blob, nombre);
        aviso(bd, "Descargada");
      });
    };
    caja.appendChild(bd);

    return caja;
  }

  function nombreArchivo(entrada, tipo) {
    var base = (entrada || "formula").toString().toLowerCase()
      .replace(/[áàäâã]/g, "a").replace(/[éèëê]/g, "e")
      .replace(/[íìïî]/g, "i").replace(/[óòöô]/g, "o")
      .replace(/[úùüû]/g, "u").replace(/ñ/g, "n")
      .replace(/[^a-z0-9]+/g, "-").replace(/^-+|-+$/g, "");
    return base + "-" + tipo + ".png";
  }

  /* ------------------------------------------------ compuestos no listados */

  function pintarNoEncontrado(consulta) {
    elDetalle.hidden = true;
    elVacio.hidden = true;
    elNoEncontrado.hidden = false;
    elNoEncontrado.className = "tarjeta aviso-caja";
    elNoEncontrado.innerHTML = "";

    var h = document.createElement("h3");
    h.textContent = "«" + consulta + "» no está en los apuntes";
    elNoEncontrado.appendChild(h);

    var parecidos = buscar(consulta).slice(0, 6);
    if (parecidos.length) {
      var p = document.createElement("p");
      p.textContent = "Quizá buscabas:";
      elNoEncontrado.appendChild(p);
      var ul = document.createElement("ul");
      ul.style.margin = "4px 0 12px";
      ul.style.paddingLeft = "20px";
      parecidos.forEach(function (r) {
        var li = document.createElement("li");
        var a = document.createElement("a");
        a.href = "#";
        a.textContent = r.e.nombre;
        a.style.color = "var(--acento)";
        a.onclick = function (ev) { ev.preventDefault(); mostrar(r.e); };
        li.appendChild(a);
        ul.appendChild(li);
      });
      elNoEncontrado.appendChild(ul);
    }

    var p2 = document.createElement("p");
    p2.textContent = "Puedo intentar construir la fórmula de todas formas. " +
      "El resultado es automático y puede contener errores: revísalo antes de fiarte.";
    elNoEncontrado.appendChild(p2);

    var bts = document.createElement("div");
    bts.className = "botones";
    var b = document.createElement("button");
    b.type = "button";
    b.className = "boton principal";
    b.textContent = "Intentar construir «" + consulta + "»";
    b.onclick = function () { construir(consulta); };
    bts.appendChild(b);
    elNoEncontrado.appendChild(bts);
  }

  function construir(nombre) {
    var r = window.Parser.construir(nombre);
    if (!r.ok) {
      var e = document.createElement("div");
      e.className = "error";
      e.textContent = "No he podido interpretarlo: " + r.error +
        " Prueba a escribirlo de otra forma o añádelo al archivo de datos.";
      elNoEncontrado.appendChild(e);
      return;
    }
    if (!r.smiles) {
      var e2 = document.createElement("div");
      e2.className = "error";
      e2.textContent = "No he podido construir la estructura.";
      elNoEncontrado.appendChild(e2);
      return;
    }

    var caja = document.createElement("div");
    caja.className = "tarjeta generado";
    caja.style.marginTop = "14px";
    caja.innerHTML = "";
    var cab = document.createElement("div");
    cab.className = "detalle-cabecera";
    var izq = document.createElement("div");
    var h = document.createElement("h2");
    h.className = "detalle-nombre";
    h.textContent = r.nombre || nombre;
    izq.appendChild(h);
    var meta = document.createElement("div");
    meta.className = "detalle-meta";
    var sello = document.createElement("span");
    sello.className = "sello";
    sello.textContent = "generado automáticamente · no está en los apuntes";
    meta.appendChild(sello);
    var fmol = document.createElement("span");
    fmol.className = "formula-molecular";
    meta.appendChild(fmol);
    izq.appendChild(meta);
    cab.appendChild(izq);
    caja.appendChild(cab);

    if (r.avisos && r.avisos.length) {
      var av = document.createElement("div");
      av.className = "notas";
      av.textContent = r.avisos.join(" ");
      caja.appendChild(av);
    }

    var vistas = document.createElement("div");
    vistas.className = "vistas";
    var v1 = document.createElement("div");
    v1.className = "vista";
    v1.innerHTML = '<div class="titulo">Estructura</div><div class="zona"></div>';
    vistas.appendChild(v1);
    caja.appendChild(vistas);
    elNoEncontrado.appendChild(caja);

    var contenedor = v1.querySelector(".zona");
    if (r.condensada) {
      var v2 = document.createElement("div");
      v2.className = "vista";
      var tt = document.createElement("div");
      tt.className = "titulo";
      tt.textContent = "Fórmula semidesarrollada";
      v2.appendChild(tt);
      var pre = document.createElement("div");
      pre.style.fontFamily = '"Cambria", "Georgia", serif';
      pre.style.fontSize = "19px";
      pre.style.padding = "20px 6px";
      pre.textContent = conSubindices(r.condensada);
      v2.appendChild(pre);
      vistas.appendChild(v2);
    }

    try {
      SmilesDrawer.parse(r.smiles, function (tree) {
        var svg = document.createElementNS("http://www.w3.org/2000/svg", "svg");
        svg.setAttribute("width", "440");
        svg.setAttribute("height", "320");
        svg.setAttribute("viewBox", "0 0 440 320");
        svg.style.width = "100%";
        svg.style.height = "auto";
        svg.style.maxHeight = "260px";
        contenedor.appendChild(svg);
        var d = new SmilesDrawer.SvgDrawer({
          width: 440, height: 320, bondThickness: 1.6,
          padding: 14, terminalCarbons: false, compactDrawing: false
        });
        d.draw(tree, svg, "oldschool", false);
        contenedor.appendChild(accionesImagen(svg, nombreArchivo(r.nombre, "estructura")));
        try {
          fmol.textContent = conSubindices(d.getMolecularFormula());
        } catch (err) { /* sin fórmula molecular */ }
      }, function () {
        contenedor.innerHTML = '<div class="titulo">Estructura</div>' +
          '<p class="error">No he podido dibujar la estructura.</p>';
      });
    } catch (err) {
      contenedor.innerHTML = '<div class="titulo">Estructura</div>' +
        '<p class="error">No he podido dibujar la estructura.</p>';
    }
  }

  /* --------------------------------------------------------- sugerencias */

  function pintarSugerencias(consulta) {
    var res = buscar(consulta).slice(0, 12);
    sugerencias = res;
    sugerenciaActiva = -1;
    if (!res.length || clave(consulta).length < 2) {
      elSug.hidden = true;
      elSug.innerHTML = "";
      return;
    }
    elSug.innerHTML = "";
    res.forEach(function (r, i) {
      var li = document.createElement("li");
      var sp = document.createElement("span");
      sp.textContent = r.e.nombre;
      li.appendChild(sp);
      var fam = document.createElement("span");
      fam.className = "fam";
      fam.textContent = familiasDe(r.e).join(" · ");
      li.appendChild(fam);
      li.onclick = function () {
        elQ.value = r.e.nombre;
        elSug.hidden = true;
        mostrar(r.e);
      };
      li.onmouseenter = function () { marcarSugerencia(i); };
      elSug.appendChild(li);
    });
    elSug.hidden = false;
  }

  function marcarSugerencia(i) {
    sugerenciaActiva = i;
    var hijos = elSug.children;
    for (var k = 0; k < hijos.length; k++) {
      hijos[k].className = (k === i) ? "activo" : "";
    }
  }

  /* ---------------------------------------------------------------- eventos */

  var temporizador = null;
  elQ.addEventListener("input", function () {
    var v = elQ.value.trim();
    elLimpiar.className = v ? "visible" : "";
    clearTimeout(temporizador);
    temporizador = setTimeout(function () {
      if (!v) {
        elSug.hidden = true;
        elDetalle.hidden = true;
        elNoEncontrado.hidden = true;
        elVacio.hidden = false;
        return;
      }
      var res = buscar(v);
      pintarSugerencias(v);
      if (!res.length) {
        pintarNoEncontrado(v);
      } else {
        elNoEncontrado.hidden = true;
        if (!elDetalle.hidden) return; /* no cambiamos lo que está mirando */
      }
    }, 90);
  });

  elQ.addEventListener("keydown", function (ev) {
    if (ev.key === "ArrowDown" && !elSug.hidden && sugerencias.length) {
      ev.preventDefault();
      marcarSugerencia(Math.min(sugerenciaActiva + 1, sugerencias.length - 1));
    } else if (ev.key === "ArrowUp" && !elSug.hidden && sugerencias.length) {
      ev.preventDefault();
      marcarSugerencia(Math.max(sugerenciaActiva - 1, 0));
    } else if (ev.key === "Enter") {
      var elegido = sugerencias[sugerenciaActiva >= 0 ? sugerenciaActiva : 0];
      if (elegido) {
        elQ.value = elegido.e.nombre;
        elSug.hidden = true;
        mostrar(elegido.e);
      }
    } else if (ev.key === "Escape") {
      elSug.hidden = true;
    }
  });

  elQ.addEventListener("focus", function () {
    if (elQ.value.trim()) pintarSugerencias(elQ.value.trim());
  });

  document.addEventListener("click", function (ev) {
    if (!elSug.contains(ev.target) && ev.target !== elQ) elSug.hidden = true;
  });

  elLimpiar.onclick = function () {
    elQ.value = "";
    elLimpiar.className = "";
    elSug.hidden = true;
    elDetalle.hidden = true;
    elNoEncontrado.hidden = true;
    elVacio.hidden = false;
    elQ.focus();
  };

  /* ------------------------------------------------------------ arranque */

  function aplicarHash() {
    var m = (location.hash || "").match(/^#(q|id)=(.+)$/);
    if (!m) return;
    var valor = decodeURIComponent(m[2].replace(/\+/g, " "));
    if (m[1] === "id") {
      for (var i = 0; i < entradas.length; i++) {
        if (entradas[i].id === valor) {
          elQ.value = entradas[i].nombre;
          elLimpiar.className = "visible";
          mostrar(entradas[i]);
          return;
        }
      }
    }
    elQ.value = valor;
    elLimpiar.className = "visible";
    var res = buscar(valor);
    if (res.length) {
      pintarSugerencias(valor);
      mostrar(res[0].e);
    } else {
      pintarNoEncontrado(valor);
    }
  }

  pintarFiltros();
  pintarLista();
  aplicarHash();
  elQ.focus();
})();
