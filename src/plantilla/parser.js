/* ==========================================================================
   Parser de nombres de química orgánica (español) -> SMILES
   ==========================================================================
   Interpreta nombres sistemáticos del tipo:

     3-metil-1-buteno        2,4-dimetilpentan-3-ona
     ácido 3-oxobutanoico    4-hidroxi-4-metil-2-pentanona
     1-metilciclopenteno     7-bromo-1,3,5-cicloheptatrieno
     1,2-diclorobenceno      2-etoxietanol

   Devuelve  { ok, smiles, condensada, error }.
   Si no entiende el nombre, ok = false y explica por qué. Nunca inventa.

   Funciona igual en el navegador (window.Parser) y en Node (module.exports).
   ========================================================================== */
(function (raiz) {
  "use strict";

  /* ------------------------------------------------------------------ */
  /*  Utilidades                                                         */
  /* ------------------------------------------------------------------ */

  function normalizar(texto) {
    return String(texto || "")
      .toLowerCase()
      .replace(/[áàäâã]/g, "a")
      .replace(/[éèëê]/g, "e")
      .replace(/[íìïî]/g, "i")
      .replace(/[óòöô]/g, "o")
      .replace(/[úùüû]/g, "u")
      .replace(/ñ/g, "n")
      // Localizadores griegos de los apuntes: α = 2, β = 3...
      .replace(/α/g, "2-")
      .replace(/β/g, "3-")
      .replace(/γ/g, "4-")
      .replace(/δ/g, "5-")
      .replace(/ε/g, "6-")
      .replace(/[·•]/g, "-")
      .replace(/\(|\)|;/g, "")
      .replace(/\s*-\s*/g, "-")
      .replace(/-+/g, "-")
      .replace(/\s+/g, " ")
      .trim();
  }

  function separar(texto) {
    return texto.split(/[-\s]+/)
      .map(function (p) { return p.trim(); })
      .filter(function (p) { return p !== ""; });
  }

  function esNumero(parte) {
    return /^\d+(,\d+)*$/.test(parte);
  }

  function numeros(parte) {
    return parte.split(",").map(Number);
  }

  /* ------------------------------------------------------------------ */
  /*  Tablas                                                             */
  /* ------------------------------------------------------------------ */

  var RAICES = {
    met: 1, et: 2, prop: 3, but: 4, pent: 5, hex: 6, hept: 7,
    oct: 8, non: 9, dec: 10, undec: 11, dodec: 12
  };
  var RAICES_ORDEN = Object.keys(RAICES).sort(function (a, b) {
    return b.length - a.length;
  });

  // Sufijos funcionales (los más largos primero)
  var SUFIJOS = [
    { clave: "amida", tipo: "cadena", fragmento: "(=O)(N)" },
    { clave: "amina", tipo: "grupo", fragmento: "(N)" },
    { clave: "oico", tipo: "cadena", fragmento: "(=O)(O)" },
    { clave: "oilo", tipo: "cadena", fragmento: "(=O)(Cl)" },
    { clave: "ona", tipo: "grupo", fragmento: "(=O)" },
    { clave: "ol", tipo: "grupo", fragmento: "(O)" },
    { clave: "al", tipo: "cadena", fragmento: "(=O)" }
  ];

  var MULTIPLICADORES = { di: 2, tri: 3, tetra: 4, penta: 5, hexa: 6 };

  var SUSTITUYENTES = {
    metil: "C",
    etil: "CC",
    propil: "CCC",
    isopropil: "C(C)C",
    "1-metiletil": "C(C)C",
    butil: "CCCC",
    isobutil: "CC(C)C",
    isopentil: "CCC(C)C",
    secbutil: "C(C)CC",
    "sec-butil": "C(C)CC",
    "tert-butil": "C(C)(C)C",
    terbutil: "C(C)(C)C",
    "terc-butil": "C(C)(C)C",
    pentil: "CCCCC",
    fenil: "c1ccccc1",
    bencil: "Cc1ccccc1",
    vinil: "C=C",
    etenil: "C=C",
    alilo: "CC=C",
    "2-propenil": "CC=C",
    "1-propenil": "C=CC",
    propenil: "CC=C",
    etinil: "C#C",
    etilnil: "C#C",
    ciclopropil: "C1CC1",
    ciclobutil: "C1CCC1",
    ciclopentil: "C1CCCC1",
    ciclohexil: "C1CCCCC1",
    fluoro: "F",
    cloro: "Cl",
    bromo: "Br",
    yodo: "I",
    hidroxi: "O",
    formil: "C=O",
    oxo: "=O",
    nitro: "[N+](=O)[O-]",
    amino: "N",
    metoxi: "OC",
    etoxi: "OCC",
    propoxi: "OCCC",
    clorometil: "CCl",
    clorometoxi: "OCCl",
    cloroperoxi: "OOC=O",
    peroxi: "OO",
    ciano: "C#N"
  };

  // Nombres triviales que no se pueden deducir de las reglas
  var TRIVIALES = {
    "fenol": "Oc1ccccc1",
    "anilina": "Nc1ccccc1",
    "tolueno": "Cc1ccccc1",
    "anisol": "COc1ccccc1",
    "estireno": "C=Cc1ccccc1",
    "xileno": "Cc1cccc(C)c1",
    "m-xileno": "Cc1cccc(C)c1",
    "meta-xileno": "Cc1cccc(C)c1",
    "o-xileno": "Cc1ccccc1C",
    "orto-xileno": "Cc1ccccc1C",
    "p-xileno": "Cc1ccc(C)cc1",
    "para-xileno": "Cc1ccc(C)cc1",
    "mesitileno": "Cc1cc(C)cc(C)c1",
    "acetofenona": "CC(=O)c1ccccc1",
    "benzaldehido": "O=Cc1ccccc1",
    "acido benzoico": "OC(=O)c1ccccc1",
    "acido formico": "OC=O",
    "acido acetico": "CC(=O)O",
    "eter dietilico": "CCOCC",
    "eter etilico": "CCOCC",
    "dietil eter": "CCOCC",
    "eter difenilico": "c1ccc(Oc2ccccc2)cc1",
    "fenoxibenceno": "c1ccc(Oc2ccccc2)cc1",
    "metano": "C",
    "etano": "CC",
    "propano": "CCC",
    "butano": "CCCC",
    "pentano": "CCCCC",
    "hexano": "CCCCCC",
    "ciclohexano": "C1CCCCC1",
    "ciclopentano": "C1CCCC1",
    "ciclobutano": "C1CCC1",
    "ciclopropano": "C1CC1",
    "benceno": "c1ccccc1",
    "metil fenil cetona": "CC(=O)c1ccccc1",
    "fenilmetilcetona": "CC(=O)c1ccccc1",
    "metil etil cetona": "CCC(C)=O",
    "etil metil cetona": "CCC(C)=O",
    "diisopropil cetona": "CC(C)C(=O)C(C)C",
    "metil neopentil cetona": "CC(=O)CC(C)(C)C",
    "propiofenona": "CCC(=O)c1ccccc1",
    "fenil etil cetona": "CCC(=O)c1ccccc1",
    "diacetona alcohol": "CC(=O)CC(C)(C)O",
    "aldehido benzoico": "O=Cc1ccccc1",
    "acido acetoacetico": "CC(=O)CC(=O)O",
    "acido 3-oxobutirico": "CC(=O)CC(=O)O",
    "aldol": "CC(O)CC=O",
    "mcpba": "OOC(=O)c1cccc(Cl)c1",
    "fenilmetilacetileno": "CC#CCc1ccccc1",
    "eter fenilico": "c1ccc(Oc2ccccc2)cc1",
    "difenil eter": "c1ccc(Oc2ccccc2)cc1",
    "alcohol feniletilico": "OCCc1ccccc1",
    "eter etil metilico": "COCC",
    "etil metil eter": "COCC",
    "clorometil metil eter": "ClCOC",
    "tetrahidrofurano": "C1CCOC1",
    "thf": "C1CCOC1",
    "oxolano": "C1CCOC1",
    "eter monoetilico del etilenglicol": "CCOCCO",
    "ciclohexano con aldehido": "O=CC1CCCCC1",
    "etileno": "C=C",
    "propileno": "C=CC",
    "acetileno": "C#C",
    "metilacetileno": "CC#C",
    "dimetilacetileno": "CC#CC",
    "fenilacetileno": "C#Cc1ccccc1",
    "feniletileno": "C=Cc1ccccc1",
    "fenilamina": "Nc1ccccc1",
    "fenil metil eter": "COc1ccccc1",
    "eter fenilmetilico": "COc1ccccc1",
    "meta-metiltolueno": "Cc1cccc(C)c1",
    "trinitrobenceno": "O=[N+]([O-])c1cc([N+](=O)[O-])cc([N+](=O)[O-])c1",
    "agua": "O",
    "metanol": "CO",
    "etanol": "CCO",
    "piperidina": "C1CCNCC1",
    "difenilamina": "N(c1ccccc1)c1ccccc1",
    "n-fenilanilina": "N(c1ccccc1)c1ccccc1",
    "n-etilanilina": "CCNc1ccccc1",
    "n,n-dietilanilina": "CCN(CC)c1ccccc1",
    "4-metilanilina": "Cc1ccc(N)cc1",
    "p-metilanilina": "Cc1ccc(N)cc1",
    "p-toluidina": "Cc1ccc(N)cc1",
    "benzoato de fenilo": "O=C(Oc1ccccc1)c1ccccc1",
    "acetato de etilo": "CCOC(C)=O",
    "acetato de metilo": "COC(C)=O",
    "formiato de isopropilo": "CC(C)OC=O",
    "fenilacetato de metilo": "COC(=O)Cc1ccccc1",
    "acido propionico": "CCC(=O)O",
    "acido oxalico": "O=C(O)C(=O)O",
    "acido glutarico": "O=C(O)CCCC(=O)O",
    "acido ftalico": "O=C(O)c1ccccc1C(=O)O",
    "acido isoftalico": "O=C(O)c1cccc(C(=O)O)c1",
    "acido m-ftalico": "O=C(O)c1cccc(C(=O)O)c1",
    "acido benceno-1,2-dicarboxilico": "O=C(O)c1ccccc1C(=O)O",
    "acido benceno-1,3-dicarboxilico": "O=C(O)c1cccc(C(=O)O)c1",
    "acido benzoato": "[O-]C(=O)c1ccccc1",
    "bencenoato": "[O-]C(=O)c1ccccc1",
    "benzoato": "[O-]C(=O)c1ccccc1",
    "ion benzoato": "[O-]C(=O)c1ccccc1",
    "anion benzoato": "[O-]C(=O)c1ccccc1",
    "gaba": "NCCCC(=O)O",
    "paba": "Nc1ccc(C(=O)O)cc1",
    "dmf": "CN(C)C=O"
  };

  /* ------------------------------------------------------------------ */
  /*  Búsqueda de la raíz principal                                      */
  /* ------------------------------------------------------------------ */

  var CONTINUACION = /^(an|en|in|a|o|ol|ona|al|oico|amida|amina|oato|ano|eno|ino)/;

  function es_continuacion_valida(resto) {
    // Quita separadores y localizadores iniciales: "-3-ino" -> "ino"
    var limpio = resto.replace(/^-/, "").replace(/^\d+(,\d+)*/, "").replace(/^-/, "");
    if (!CONTINUACION.test(limpio)) return false;
    // "metil", "metoxi", "etinil", "propenil"... son sustituyentes
    if (/^(il|oxi|ilo|enil|anil|inil)/.test(limpio)) return false;
    return true;
  }

  function buscar_raiz(texto) {
    var especiales = [
      { clave: "bencen", n: 6 },
      { clave: "benz", n: 6 },
      { clave: "fen", n: 6 }
    ];
    for (var i = 0; i < texto.length; i++) {
      for (var j = 0; j < especiales.length; j++) {
        var e = especiales[j];
        if (texto.substr(i, e.clave.length) === e.clave) {
          var resto = texto.substr(i + e.clave.length);
          if (resto === "" || es_continuacion_valida(resto)) {
            return {
              inicio: i, fin: i + e.clave.length, n: e.n,
              aromatico: true, anillo: true, clave: e.clave
            };
          }
        }
      }
      for (var k = 0; k < RAICES_ORDEN.length; k++) {
        var clave = RAICES_ORDEN[k];
        if (texto.substr(i, clave.length) === clave) {
          var resto2 = texto.substr(i + clave.length);
          if (es_continuacion_valida(resto2)) {
            var antes = texto.substr(0, i);
            var ciclo = /ciclo$/.test(antes);
            return {
              inicio: ciclo ? i - 5 : i,
              fin: i + clave.length,
              n: RAICES[clave],
              aromatico: false,
              anillo: ciclo,
              clave: clave
            };
          }
        }
      }
    }
    return null;
  }

  /* ------------------------------------------------------------------ */
  /*  Interpretación del bloque derecho (enlaces y sufijo funcional)     */
  /* ------------------------------------------------------------------ */

  function extraer_localizador(texto) {
    var m = texto.match(/(\d+(?:,\d+)*)-?$/);
    if (!m) return { locantes: [], resto: texto };
    return {
      locantes: numeros(m[1]),
      resto: texto.substr(0, texto.length - m[0].length)
    };
  }

  function parsear_sufijo(texto) {
    var resto = texto;
    var res = { enlaces: [], sufijo: null };

    // 0) Diácidos: «etanodioico», «pentanodioico» (dos -COOH)
    if (/dioico$/.test(resto)) {
      var antesDi = resto.substr(0, resto.length - "dioico".length);
      var locDi = extraer_localizador(antesDi);
      res.sufijo = {
        tipo: "dioico",
        fragmento: "(=O)(O)",
        locante: null,
        nombre: "dioico"
      };
      resto = locDi.resto;
    }

    // 1) Sufijo funcional al final
    for (var s = 0; !res.sufijo && s < SUFIJOS.length; s++) {
      var suf = SUFIJOS[s];
      var idx = resto.lastIndexOf(suf.clave);
      if (idx < 0) continue;
      var despues = resto.substr(idx + suf.clave.length);
      if (despues !== "" && despues !== "o" && despues !== "a") continue;
      var antes = resto.substr(0, idx);
      var loc = extraer_localizador(antes);
      res.sufijo = {
        tipo: suf.tipo,
        fragmento: suf.fragmento,
        locante: loc.locantes.length ? loc.locantes[0] : null,
        nombre: suf.clave
      };
      resto = loc.resto;
      break;
    }

    // 2) Enlaces múltiples: "adieno", "-2-en", "ano", "-3-ino"
    var re = /(di|tri|tetra|penta|hexa)?(en|in)/g;
    var m;
    while ((m = re.exec(resto))) {
      var mult = MULTIPLICADORES[m[1]] || 1;
      var antesEnlace = resto.substr(0, m.index);
      var locEnlace = extraer_localizador(antesEnlace);
      var locs = locEnlace.locantes;
      for (var i = 0; i < mult; i++) {
        res.enlaces.push({
          tipo: m[2],
          locante: locs[i] !== undefined ? locs[i] : null
        });
      }
      resto = locEnlace.resto + resto.substr(m.index + m[0].length);
      re.lastIndex = 0;
    }

    // 3) Lo que queda solo puede ser "an", "a", "o" o separadores
    var limpio = resto.replace(/[-\s]/g, "").replace(/^an/, "").replace(/[ao]/g, "");
    return { resultado: res, limpio: limpio };
  }

  /* ------------------------------------------------------------------ */
  /*  Interpretación del bloque izquierdo (sustituyentes)                */
  /* ------------------------------------------------------------------ */

  // Trocea un bloque sin separadores en sustituyentes:
  // «dietilmetil» -> di-etil + metil, «ciclohexildimetil» -> ciclohexil + di-metil.
  var SUSTITUYENTES_ORDEN = Object.keys(SUSTITUYENTES).sort(function (a, b) {
    return b.length - a.length;
  });

  function trocear_sustituyentes(texto) {
    var trozos = [];
    var resto = texto;
    while (resto !== "") {
      var mult = 1;
      var m = resto.match(/^(di|tri|tetra|penta|hexa)/);
      if (m) {
        mult = MULTIPLICADORES[m[1]];
        resto = resto.substr(m[1].length);
      }
      var encontrado = null;
      for (var k = 0; k < SUSTITUYENTES_ORDEN.length; k++) {
        var clave = SUSTITUYENTES_ORDEN[k];
        if (resto.substr(0, clave.length) === clave) { encontrado = clave; break; }
      }
      if (!encontrado) return null;
      trozos.push({ nombre: encontrado, mult: mult });
      resto = resto.substr(encontrado.length);
    }
    return trozos;
  }

  function parsear_bloque_izquierdo(texto) {
    var partes = separar(texto);
    var sustituyentes = [];
    var n_sustituyentes = [];
    var pendientes = [];
    var nPendientes = 0;
    var i = 0;
    while (i < partes.length) {
      var parte = partes[i];
      if (esNumero(parte)) {
        pendientes = numeros(parte);
        i++;
        continue;
      }
      // Localizadores del nitrógeno: «N-», «N,N-», «N,2-»
      if (/^n(,(\d+|n))*$/.test(parte)) {
        var trozosLoc = parte.split(",");
        nPendientes = 0;
        var numerosLoc = [];
        for (var q = 0; q < trozosLoc.length; q++) {
          if (trozosLoc[q] === "n") nPendientes++;
          else numerosLoc.push(Number(trozosLoc[q]));
        }
        pendientes = numerosLoc.concat(pendientes);
        i++;
        continue;
      }
      if (/^(tert|terc|sec|iso)$/.test(parte) && i + 1 < partes.length) {
        parte = parte + "-" + partes[i + 1];
        i++;
      }
      var lecturas = trocear_sustituyentes(parte);
      if (!lecturas) {
        return { error: "No reconozco el sustituyente «" + parte + "»." };
      }
      var locantes = pendientes.slice();
      pendientes = [];
      var usados = 0;
      for (var t = 0; t < lecturas.length; t++) {
        for (var k = 0; k < lecturas[t].mult; k++) {
          var alNitrogeno = nPendientes > 0;
          if (alNitrogeno) nPendientes--;
          var locante = null;
          if (!alNitrogeno) {
            locante = locantes[usados] !== undefined ? locantes[usados]
                     : (locantes[0] !== undefined ? locantes[0] : null);
            usados++;
          }
          var su = {
            nombre: lecturas[t].nombre,
            fragmento: SUSTITUYENTES[lecturas[t].nombre],
            locante: locante
          };
          if (alNitrogeno) n_sustituyentes.push(su);
          else sustituyentes.push(su);
        }
      }
      i++;
    }
    return {
      sustituyentes: sustituyentes,
      n_sustituyentes: n_sustituyentes,
      locantes_sueltos: pendientes
    };
  }

  /* ------------------------------------------------------------------ */
  /*  Construcción del esqueleto                                         */
  /* ------------------------------------------------------------------ */

  function renumerarCierres(fragmento, contador) {
    if (!/\d/.test(fragmento)) return { fragmento: fragmento, contador: contador };
    var mapa = {};
    var nuevo = fragmento.replace(/\d/g, function (d) {
      if (!mapa[d]) {
        mapa[d] = contador > 9 ? "%" + contador : String(contador);
        contador++;
      }
      return mapa[d];
    });
    return { fragmento: nuevo, contador: contador };
  }

  function construirTokens(raiz, enlaces) {
    var n = raiz.n;
    var tokens = [];
    var i;
    if (raiz.aromatico) {
      for (i = 0; i < 6; i++) tokens.push(i === 0 ? "c1" : (i === 5 ? "c1" : "c"));
    } else if (raiz.anillo) {
      for (i = 0; i < n; i++) tokens.push(i === 0 ? "C1" : (i === n - 1 ? "C1" : "C"));
    } else {
      for (i = 0; i < n; i++) tokens.push("C");
    }
    // Enlaces múltiples: el localizador k indica el enlace entre los átomos
    // k-1 y k (1-based); el símbolo se escribe delante del segundo átomo.
    for (i = 0; i < enlaces.length; i++) {
      var e = enlaces[i];
      if (e.locante === null) continue;
      var idx = e.locante % n; // locante n -> enlace de cierre
      tokens[idx] = (e.tipo === "en" ? "=" : "#") + tokens[idx];
    }
    return tokens;
  }

  /* ------------------------------------------------------------------ */
  /*  Fórmula semidesarrollada (texto) para compuestos generados         */
  /* ------------------------------------------------------------------ */

  var NOMBRES_RAMA = {
    "C": "CH3", "CC": "CH2CH3", "CCC": "CH2CH2CH3", "O": "OH",
    "=O": "=O", "C=O": "CHO", "(=O)(O)": "COOH", "(=O)(OO)": "COOOH",
    "N": "NH2", "OC": "OCH3", "OCC": "OCH2CH3", "Cl": "Cl", "Br": "Br",
    "F": "F", "I": "I", "c1ccccc1": "Ph"
  };

  function nombreRama(fragmento) {
    return NOMBRES_RAMA[fragmento] || fragmento;
  }

  function condensadaDesdeTokens(tokens) {
    var n = tokens.length;
    var salida = [];
    for (var i = 0; i < n; i++) {
      var t = tokens[i];
      var prefijo = (t.match(/^[=#]/) || [""])[0];
      var ramas = [];
      var re = /\(([^()]*)\)/g;
      var m;
      while ((m = re.exec(t))) ramas.push(m[1]);
      var base = t.replace(/^[=#]/, "").replace(/\([^()]*\)/g, "").replace(/1/g, "");

      // enlaces de este átomo
      var enlaces = 0;
      if (i > 0) enlaces += prefijo === "=" ? 2 : prefijo === "#" ? 3 : 1;
      if (i < n - 1) {
        var sig = (tokens[i + 1].match(/^[=#]/) || [""])[0];
        enlaces += sig === "=" ? 2 : sig === "#" ? 3 : 1;
      }
      for (var r = 0; r < ramas.length; r++) {
        enlaces += ramas[r].charAt(0) === "=" ? 2 : 1;
      }
      var hidrogenos = Math.max(0, 4 - enlaces);
      var grupo = base + (hidrogenos === 0 ? "" : hidrogenos === 1 ? "H" : "H" + hidrogenos);

      // casos especiales legibles
      var clave = ramas.join("");
      if (ramas.length === 2 && ramas.indexOf("=O") >= 0 && ramas.indexOf("O") >= 0) {
        grupo = "COOH";
      } else if (ramas.length === 1 && ramas[0] === "=O" && hidrogenos === 1) {
        grupo = "CHO";
      } else if (ramas.length) {
        grupo += "(" + ramas.map(nombreRama).join(",") + ")";
      }
      var separador = "";
      if (i > 0) {
        var pre = (t.match(/^[=#]/) || [""])[0];
        separador = pre === "=" ? "=" : pre === "#" ? "#" : "-";
      }
      salida.push(separador + grupo);
    }
    return salida.join("");
  }

  /* ------------------------------------------------------------------ */
  /*  Aminas nombradas por sus sustituyentes                             */
  /* ------------------------------------------------------------------ */

  function amina_desde_sustituyentes(texto) {
    // «etilamina», «bencilamina», «ciclohexilamina», «tert-butilamina»,
    // «dietilmetilamina», «ciclohexildimetilamina», «difenilamina»...
    if (!/amina$/.test(texto)) return null;
    var bloque = texto.replace(/amina$/, "").replace(/-$/, "");
    if (!bloque) return null;
    var pi = parsear_bloque_izquierdo(bloque);
    if (pi.error || !pi.sustituyentes.length || pi.locantes_sueltos.length) return null;
    var contador = 2;
    var smiles = "N";
    for (var i = 0; i < pi.sustituyentes.length; i++) {
      var f = renumerarCierres(pi.sustituyentes[i].fragmento, contador);
      contador = f.contador;
      smiles += "(" + f.fragmento + ")";
    }
    return { ok: true, smiles: smiles, condensada: "", avisos: [], nombre: texto };
  }

  /* ------------------------------------------------------------------ */
  /*  Función principal                                                  */
  /* ------------------------------------------------------------------ */

  function construir(nombreOriginal, opciones) {
    var texto = normalizar(nombreOriginal);
    if (!texto) return { ok: false, error: "Escribe un nombre." };
    opciones = opciones || {};

    if (TRIVIALES[texto]) {
      return { ok: true, smiles: TRIVIALES[texto], condensada: "", trivial: true, nombre: texto };
    }

    // Ésteres: «etanoato de etilo» = ácido etanoico + etilo
    var mEster = texto.match(/^(.+?)ato de (.+)$/);
    if (mEster) {
      var alquilo = SUSTITUYENTES[mEster[2].replace(/o$/, "")];
      if (alquilo) {
        var rEster = construir("acido " + mEster[1] + "ico", { ester: alquilo });
        if (rEster.ok) {
          rEster.nombre = texto;
          rEster.condensada = "";
          return rEster;
        }
      }
    }

    var esAcido = /^acido\s+/.test(texto);
    if (esAcido) texto = texto.replace(/^acido\s+/, "");

    var posicion = 0;
    // Prefijo de posición (orto/meta/para). Ojo: «metanoico», «metanoato»...
    // empiezan por «meta» pero no son prefijos de posición.
    var rePosicion = /(^|[- ])(orto|meta|para)(?=(?:-|[^n]|n[^ao]))|(^|[- ])(o|m|p)(?=[- ])/;
    var mPos = texto.match(rePosicion);
    if (mPos) {
      var palabra = mPos[2] || mPos[4];
      posicion = (palabra === "orto" || palabra === "o") ? 2
               : (palabra === "meta" || palabra === "m") ? 3 : 4;
      texto = texto.replace(rePosicion,
                            function (todo, a, b, c, d) { return (a === "-" || c === "-") ? "-" : ""; });
    }

    var carbaldehido = /carb(ox)?aldehido/.test(texto);
    if (carbaldehido) {
      texto = texto.replace(/-?carb(ox)?aldehido/, "");
    }
    var carboxamida = /carbox?amida/.test(texto);
    if (carboxamida) {
      texto = texto.replace(/-?carbox?amida/, "");
    }
    var carboxilico = /carboxilico/.test(texto);
    if (carboxilico) {
      texto = texto.replace(/-?carboxilico/, "");
    }

    var raiz = buscar_raiz(texto);
    if (!raiz) {
      // «etilamina», «bencilamina», «dietilmetilamina»: el nitrógeno cuelga
      // de los sustituyentes que se nombran en el bloque izquierdo.
      var aminaDirecta = amina_desde_sustituyentes(texto);
      if (aminaDirecta) return aminaDirecta;
      return { ok: false, error: "No reconozco ninguna cadena principal en «" + nombreOriginal + "»." };
    }

    var textoIzq = texto.substr(0, raiz.inicio).replace(/-$/, "");
    var textoDer = texto.substr(raiz.fin);

    var avisos = [];
    var contadorCierres = 2;

    // Peroxiácido: "peroxibenzoico" -> el grupo ácido lleva un O de más
    var peroxi = /peroxi/.test(textoIzq);
    if (peroxi) {
      textoIzq = textoIzq.replace(/-?peroxi/, "");
      avisos.push("Es un peroxiácido: el grupo ácido lleva un oxígeno extra (COOOH).");
    }

    var ps = parsear_sufijo(textoDer);
    if (ps.limpio !== "") {
      return { ok: false, error: "No entiendo «" + ps.limpio + "» en «" + nombreOriginal + "»." };
    }
    var enlaces = ps.resultado.enlaces;
    var sufijo = ps.resultado.sufijo;

    if (esAcido && (!sufijo || sufijo.nombre !== "dioico")) {
      sufijo = {
        tipo: "cadena",
        fragmento: "(=O)(O" + (opciones.ester ? opciones.ester : "") + ")",
        locante: null,
        nombre: "oico"
      };
    }
    if (peroxi && sufijo) {
      sufijo.fragmento = "(=O)(OO)";
    }

    var pi = parsear_bloque_izquierdo(textoIzq);
    if (pi.error) return { ok: false, error: pi.error };

    // Sustituyentes del nitrógeno: «N-metil...», «N,N-dietil...»
    var nfrag = "";
    for (var n = 0; n < pi.n_sustituyentes.length; n++) {
      var fn = renumerarCierres(pi.n_sustituyentes[n].fragmento, contadorCierres);
      contadorCierres = fn.contador;
      nfrag += "(" + fn.fragmento + ")";
    }
    if (nfrag) {
      if (!sufijo || (sufijo.nombre !== "amina" && sufijo.nombre !== "amida")) {
        return { ok: false, error: "El nitrógeno sustituido solo se entiende con -amina o -amida." };
      }
      if (sufijo.nombre === "amina") {
        sufijo.fragmento = "(N" + nfrag + ")";
      } else {
        sufijo.fragmento = "(=O)(N" + nfrag + ")";
      }
    }

    // Reparto de localizadores sueltos: primero los enlaces, luego el sufijo
    var sueltos = pi.locantes_sueltos.slice();
    for (var e = 0; e < enlaces.length; e++) {
      if (enlaces[e].locante === null && sueltos.length) {
        enlaces[e].locante = sueltos.shift();
      }
    }
    if (sufijo && sufijo.locante === null && sueltos.length) {
      sufijo.locante = sueltos.shift();
    }
    // Enlaces sin localizador: 1, 3, 5... (lo habitual en dienos y trienos)
    var siguiente = 1;
    for (var e2 = 0; e2 < enlaces.length; e2++) {
      if (enlaces[e2].locante === null) {
        enlaces[e2].locante = siguiente;
        if (raiz.n >= 3 && !raiz.anillo) {
          avisos.push("He supuesto el enlace " + (enlaces[e2].tipo === "en" ? "doble" : "triple") +
                      " en el carbono " + siguiente + ".");
        }
      }
      siguiente = enlaces[e2].locante + 2;
    }

    var tokens = construirTokens(raiz, enlaces);

    // Sufijo funcional
    if (sufijo) {
      if (sufijo.tipo === "grupo" && sufijo.locante === null) {
        if (raiz.anillo) {
          sufijo.locante = 1;
        } else {
          var ocupados = {};
          for (var q = 0; q < pi.sustituyentes.length; q++) {
            if (pi.sustituyentes[q].locante) ocupados[pi.sustituyentes[q].locante] = true;
          }
          var esCetona = sufijo.nombre === "ona";
          if (esCetona ? ocupados[1] : !ocupados[1]) {
            sufijo.locante = 1;
          } else {
            var libre = 2;
            while (ocupados[libre] && libre <= tokens.length) libre++;
            sufijo.locante = libre;
          }
        }
      }
      if (sufijo.tipo === "dioico") {
        // un -COOH en cada extremo de la cadena
        tokens[0] = tokens[0] + sufijo.fragmento;
        tokens[tokens.length - 1] = tokens[tokens.length - 1] + sufijo.fragmento;
      } else if (raiz.aromatico && sufijo.tipo === "cadena") {
        // el grupo carboxilo / aldehído cuelga del carbono 1 del anillo
        tokens[0] = tokens[0] + "(C" + sufijo.fragmento + ")";
      } else {
        var loc = sufijo.locante === null ? 1 : sufijo.locante;
        if (sufijo.tipo === "cadena") loc = 1;
        if (loc < 1 || loc > tokens.length) {
          return { ok: false, error: "El localizador " + loc + " se sale de la cadena." };
        }
        tokens[loc - 1] = tokens[loc - 1] + sufijo.fragmento;
      }
    }

    // Sustituyentes
    for (var s = 0; s < pi.sustituyentes.length; s++) {
      var su = pi.sustituyentes[s];
      if (su.locante === null) {
        if (posicion && raiz.anillo) {
          var libres = 0;
          for (var z = 0; z < pi.sustituyentes.length; z++) {
            if (pi.sustituyentes[z].locante === null) libres++;
          }
          if (libres > 1 && !sufijo) {
            su.locante = s === 0 ? 1 : posicion;
          } else {
            su.locante = posicion;
          }
          avisos.push("Posición " + (su.locante === 2 ? "orto" : su.locante === 3 ? "meta" : "para") + " interpretada como carbono " + su.locante + ".");
        } else if (raiz.anillo || pi.sustituyentes.length === 1) {
          su.locante = 1;
          avisos.push("He supuesto que «" + su.nombre + "» va en el carbono 1.");
        } else {
          return { ok: false, error: "Falta el localizador del sustituyente «" + su.nombre + "»." };
        }
      }
      var idx = su.locante - 1;
      if (idx < 0 || idx >= tokens.length) {
        return { ok: false, error: "El localizador " + su.locante + " se sale de la cadena." };
      }
      var frag = renumerarCierres(su.fragmento, contadorCierres);
      contadorCierres = frag.contador;
      tokens[idx] = tokens[idx] + "(" + frag.fragmento + ")";
    }

    if (carbaldehido) {
      tokens[0] = tokens[0] + "(C=O)";
    }
    if (carboxamida) {
      tokens[0] = tokens[0] + "(C(=O)N)";
    }
    if (carboxilico && !esAcido) {
      tokens[0] = tokens[0] + "(C(=O)(O))";
    }

    var smiles = tokens.join("");
    return {
      ok: true,
      smiles: smiles,
      condensada: raiz.anillo ? "" : condensadaDesdeTokens(tokens),
      avisos: avisos,
      nombre: texto
    };
  }

  var Parser = {
    construir: construir,
    normalizar: normalizar,
    _tablas: { RAICES: RAICES, SUSTITUYENTES: SUSTITUYENTES, TRIVIALES: TRIVIALES }
  };

  if (typeof module !== "undefined" && module.exports) module.exports = Parser;
  raiz.Parser = Parser;
})(typeof window !== "undefined" ? window : globalThis);
