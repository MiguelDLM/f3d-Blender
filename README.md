# f3d-blender

Extensión de Blender para **importar archivos Autodesk Fusion 360 `.f3d`**.

Un `.f3d` es un contenedor ZIP cuya geometría está guardada como **B-rep en
formato ASM (Autodesk Shape Manager)** — no contiene malla. Esta extensión hace
ingeniería inversa del formato, parsea el B-rep y lo **tesela** a una malla de
Blender. Ver [`docs/FORMAT.md`](docs/FORMAT.md) para la documentación del
formato.

## Estado

| Componente | Fichero | Estado |
|---|---|---|
| Contenedor ZIP | `io_scene_f3d/container.py` | ✅ |
| Tokenizer ASM SAB | `io_scene_f3d/sab.py` | ✅ |
| Reconstrucción topológica B-rep | `io_scene_f3d/brep.py` | ✅ |
| Superficies analíticas (cono, extrusión, ref pool) | `io_scene_f3d/surfaces.py` | ✅ |
| Curvas B-spline (de Boor, abiertas y periódicas) | `io_scene_f3d/nurbs.py` | ✅ |
| Teselación (CDT paramétrico + loft) | `io_scene_f3d/tessellate.py` | ✅ |
| Operador de importación | `io_scene_f3d/importer.py` | ✅ |

Fidelidad contra el OBJ exportado por Fusion del mismo modelo: área total
76.668 vs 76.608 cm² (0.08 %), desviación media 10–26 µm (máx 140 µm),
malla estanca (0 aristas no-manifold).

## Pruebas (sin Blender)

```bash
python3 tests/test_parse.py
```

## Instalación en Blender (4.2+)

Empaquetar `io_scene_f3d/` como extensión e instalarla desde
*Edit > Preferences > Get Extensions > Install from Disk*.
