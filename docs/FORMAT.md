# Estructura del formato Autodesk Fusion 360 `.f3d`

Documento de ingeniería inversa del formato `.f3d` (Fusion Archive), basado en el
análisis del fichero de muestra `samples/soporte de corttinas.f3d`
(Fusion 360, escrito el 2025-01-07 con ASM 230.5.1).

> **Aviso**: `.f3d` es un formato **propietario y no documentado** de Autodesk.
> Todo lo aquí descrito procede de ingeniería inversa sobre muestras reales, no
> de documentación oficial. Es correcto para las muestras analizadas pero puede
> no cubrir todas las variantes que genera Fusion.

---

## 1. Contenedor: es un ZIP

Un `.f3d` es un archivo **ZIP** (entradas almacenadas sin comprimir,
`method=store`). Estructura de entradas observada:

```
Manifest.dat                         Metadatos del documento (texto UTF-16)
Properties.dat                       4 bytes (flags del documento)
FusionAssetName[Active]/
├── Manifest.dat                     Manifiesto del asset activo (UTF-16)
├── Breps.BlobParts/
│   └── BREP.<guid>.smb              ← GEOMETRÍA: B-rep ASM (ver §2)
├── Design1/
│   ├── BulkStream.dat               Árbol del diseño (componentes, unidades,
│   │                                 materiales, instancias) — binario propio
│   └── MetaStream.dat               Índices/metadatos del diseño
├── ProteinAssets.BlobParts/
│   └── ProteinAsset.<guid>.protein  Materiales (formato Autodesk Protein)
├── Images.BlobParts/                Texturas embebidas (si las hay)
├── DesignConfigurationTable.BlobParts/
│   └── *.dsgcfg                     Tabla de configuraciones del diseño
└── Previews/
    └── small.png                    Miniatura 256×256 RGBA
```

Los `Manifest.dat` son texto **UTF-16LE** con campos separados por bytes de
control; contienen `FusionDocType`, GUIDs del documento y tipos de asset
(`Design`, `CAM`, `SimCommon`, …).

**Dato clave:** la geometría existe **únicamente como B-rep** dentro del `.smb`.
No hay ninguna malla teselada en el archivo — Fusion la genera al vuelo para
visualizar. Por eso un importador debe parsear el B-rep y **teselarlo**.

---

## 2. Geometría: `BREP.*.smb` — ASM BinaryFile

El `.smb` empieza con la cadena mágica **`ASM BinaryFile4`** seguida de un
preámbulo fijo de 16 bytes y luego un flujo de *records* tokenizados. El
esquema es el de **ACIS "Standard ACIS Binary" (`.sab`)**: Autodesk Shape
Manager (ASM) es un derivado del kernel ACIS y comparte su serialización.

Cabecera del ejemplo:

```
ASM BinaryFile4
<16 bytes: campo de 8 bytes + dos int32>
"Autodesk Neutron"            producto
"ASM 230.5.1.65535 NT"        versión del kernel
"Tue Jan  7 12:20:18 2025"    fecha
10.0        (normalización de longitud)
1e-06       resabs  (tolerancia absoluta)
1e-10       resnor  (tolerancia normal)
0x0d "asmheader" ... 0x11(EOR)
```

### 2.1 Tokenización

El cuerpo es una secuencia de **tokens**: un byte de *tag* seguido de su carga.
Los records se delimitan **exclusivamente** por el tag de fin de record
(`0x11`). Tags identificados:

| Tag  | Significado                | Carga                                   |
|------|----------------------------|-----------------------------------------|
| `0x04` | int con signo            | 4 bytes (int32 LE)                      |
| `0x05` | float                    | 4 bytes (float32 LE) *(raro)*           |
| `0x06` | double                   | 8 bytes (IEEE-754 LE)                   |
| `0x07` | string                   | u8 longitud + bytes                     |
| `0x08` | string (uso alternativo) | u8 longitud + bytes                     |
| `0x0A` | enum / lógico            | 0 bytes                                 |
| `0x0B` | enum / lógico            | 0 bytes                                 |
| `0x0C` | **puntero** a record     | int32 = índice de record (`-1` = nulo)  |
| `0x0D` | **nombre de entidad**    | u8 longitud + nombre                    |
| `0x0E` | apertura de subtipo      | u8 longitud + nombre                    |
| `0x0F` | cierre de subtipo        | 0 bytes                                 |
| `0x10` | marcador                 | 0 bytes                                 |
| `0x11` | **fin de record (EOR)**  | 0 bytes                                 |
| `0x13` | **posición 3D**          | 3 doubles (x, y, z) = 24 bytes          |
| `0x14` | vector/dirección 3D      | 3 doubles = 24 bytes                    |
| `0x15` | uint32 (flags/contadores)| 4 bytes                                 |

**Punto crítico de la ingeniería inversa:** dentro de un mismo record puede
haber **varios** tags `0x0D`. ASM guarda las entidades de forma polimórfica
escribiendo primero el nombre de la clase base y luego el/los derivados
(p.ej. `curve` → `exact_int_cur` → `nubs`, o `surface` → `spline` → `nubs`).
Por eso `0x0D` **no** sirve como delimitador de record; sólo `0x11` lo hace.
El nombre "de tipo" útil es el **último** (la geometría concreta).

Los punteros (`0x0C`) referencian records por su **índice ordinal** en el flujo
(0-based, en orden de aparición). El parser (`io_scene_f3d/sab.py`) construye
la lista de records y resuelve estos índices.

### 2.2 Grafo topológico (B-rep)

Jerarquía clásica ACIS reconstruida a partir de los punteros:

```
body → lump → shell → face → loop → coedge → edge → vertex
                        │                       │
                        └─ surface              └─ curve
```

- **body** / **lump** / **shell**: agrupación del sólido.
- **face**: una cara delimitada; referencia una **surface** y su **loop**(s).
- **loop**: bucle cerrado de **coedge**s que recorta la cara.
- **coedge**: uso orientado de una **edge** por una cara (con `pcurve`, la curva
  en el espacio paramétrico de la superficie).
- **edge**: arista con vértices inicial/final, parámetros y una **curve**.
- **vertex** → **point**: posición 3D (tag `0x13`).

Tipos de **superficie** vistos: `plane`, `cone`, `spline` (B-spline `nubs`).
Tipos de **curva**: `straight` (recta), `ellipse`, `intcurve`/`exact_int_cur`
(curva de intersección, B-spline `nubs`).

Conteo del ejemplo `soporte de corttinas.f3d`:

| Entidad | Nº | Entidad  | Nº  | Geometría | Nº |
|---------|----|----------|-----|-----------|----|
| body    | 8  | coedge   | 648 | plane     | 64 |
| lump    | 8  | edge     | 324 | cone      | 5  |
| shell   | 8  | vertex   | 216 | spline    | 56 |
| face    | 125| loop     | 130 | straight  | 202|
| point   | 216| pcurve   | 225 | ellipse   | 8  |

### 2.3 B-spline `nubs`

Bloque de una curva/superficie B-spline no racional (`nubs`):
grado, luego los **nudos** (pares multiplicidad `int` + valor `double`) y los
**puntos de control** como ternas de `double` (x, y, z). Las `nurbs` racionales
añadirían pesos (no observadas en la muestra).

### 2.4 Atributos

Abundan los records `ATTRIB_CUSTOM` (458 en la muestra): metadatos de Fusion
adjuntos a las entidades (p.ej. `sketch_attrib_def`, IDs de tags). No aportan
geometría y el importador los ignora.

---

## 3. Estado de la ingeniería inversa

| Parte                              | Estado |
|------------------------------------|--------|
| Contenedor ZIP + manifiestos       | ✅ Documentado |
| Tokenización ASM SAB               | ✅ Parser completo (`sab.py`) |
| Grafo topológico B-rep             | ✅ Se parsea (2816 records) |
| Geometría de superficies/curvas    | 🚧 En progreso (`brep.py`) |
| Teselación de caras recortadas     | 🚧 Pendiente (`tessellate.py`) |
| `BulkStream.dat` (árbol/materiales)| 🔜 Parcial (strings) |
| Operador de importación Blender    | 🔜 Pendiente (`importer.py`) |
