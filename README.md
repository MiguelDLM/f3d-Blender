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
| Reconstrucción topológica B-rep | `io_scene_f3d/brep.py` | 🚧 |
| Teselación de caras | `io_scene_f3d/tessellate.py` | 🔜 |
| Operador de importación | `io_scene_f3d/importer.py` | 🔜 |

## Pruebas (sin Blender)

```bash
python3 tests/test_parse.py
```

## Instalación en Blender (4.2+)

Empaquetar `io_scene_f3d/` como extensión e instalarla desde
*Edit > Preferences > Get Extensions > Install from Disk*.
