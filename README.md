# EASiMap (Enzyme Active Site Map)

Mapas celestes de familias de enzimas, vistos desde el sustrato. Colócate sobre el azúcar, en el centro del enzima, y mira hacia arriba: cada estrella es un aminoácido.

EASiMap 2 es la versión nueva del trabajo de fin de máster de Miguel Ferreiro García (IQS, 2020), dirigido por el Dr. Xevi Biarnés. No se presenta como ciencia nueva. Los trabajos previos y relacionados están citados en la propia web.

## Ver la web en local

```bash
python -m http.server 8784 --bind 127.0.0.1 --directory web
```

Abre `http://127.0.0.1:8784`. No hay paso de compilación.

## Estructura

| Ruta | Contenido |
| --- | --- |
| `web/index.html` | Inicio: la metáfora, la familia GH1, las cinco familias, la idea, cómo se lee, cómo se usa, historias y método |
| `web/visor.html` | Visor de las cinco familias: mapas, filtros, selectores, Mol* y ayuda breve |
| `web/css/style.css` | Estilo |
| `web/js/i18n.js` | Todos los textos, en inglés y español |
| `web/js/core.js` | Datos, proyección y dibujo común |
| `web/js/hero.js`, `metaphor.js`, `family.js`, `skies.js`, `idea.js`, `guide.js`, `stories.js` | Piezas del inicio |
| `web/js/explorer.js`, `mol.js` | Piezas del visor |
| `web/data/gh1.json`, `gh1.extras.json` | Conjunto actual, del pipeline, y sus organismos y cifras |
| `web/data/gh5.json`, `gh16`, `gh18`, `gh20` y sus `.extras.json` | Conjuntos actuales de las otras familias |
| `web/data/families.json` | Resumen de las cinco familias para el inicio |
| `web/scripts/extras_from_pipeline.py` | Genera los extras y `families.json` |
| `web/scripts/stamp_version.py` | Pone un solo valor de versión en las dos páginas |
| `pipeline/` | Pipeline nuevo en Python |

## Datos

La web usa solo datos actuales, medidos desde el PDB con el pipeline de este repo (`source: "pipeline"`). Hay cinco familias: GH1, GH5, GH16, GH18 y GH20. CE4, la sexta del trabajo original, volverá cuando tenga un conjunto actual.

Las pruebas de réplica comparan el pipeline con la tabla de 2020. Esa tabla no se publica: las pruebas se omiten si `EASIMAP_TFM_ROOT` no apunta a ella.

## Antes de cada commit

Si cambian `web/js`, `web/css` o `web/data`, ejecuta `python web/scripts/stamp_version.py`. Los módulos, el estilo y los datos se piden con ese valor, y así un navegador no mezcla ficheros de dos subidas.

## Licencias

Fuentes de datos y software: PDB (CC0), Pfam (CC0), Mol* (MIT), D3 (ISC). Tipografías: SIL Open Font License 1.1, en `web/fonts/`.

El código de este repo se publica con licencia MIT (`LICENSE`).
