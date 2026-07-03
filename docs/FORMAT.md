# Estructura del formato Autodesk Fusion 360 `.f3d`

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

**Punto crítico:** dentro de un mismo record puede
haber **varios** tags `0x0D`. ASM guarda las entidades de forma polimórfica
escribiendo primero el nombre de la clase base y luego el/los derivados
(p.ej. `curve` → `exact_int_cur` → `nubs`, o `surface` → `spline` → `nubs`).
Por eso `0x0D` **no** sirve como delimitador de record; sólo `0x11` lo hace.
El nombre "de tipo" útil es el **último** (la geometría concreta).

Los punteros (`0x0C`) referencian records por su **índice ordinal** en el flujo
(0-based, en orden de aparición). **Ojo**: la numeración del escritor **no
cuenta** los records estructurales `delta_state` ni `Begin` (sí cuenta `End`).
En los `.smb` simples no aparecen y el ordinal coincide con la posición física;
en los `.smbh` (con diario de historial intercalado) hay que saltarlos al
resolver — validado: con esa regla el 100 % de los punteros cara→superficie
del Perchero resuelven a superficies (con numeración física, solo el 67 %).

#### `.smb` vs `.smbh`

Los `.f3d` recientes traen **dos** blobs ASM en `Breps.BlobParts/`:

- **`BREP.<guid>.smbh`** — el **modelo actual**: los cuerpos del árbol de
  diseño en sus **posiciones finales** (tras las features de mover/copiar),
  seguidos de un diario de historial (`Begin`/`delta_state`/`End`). Es el
  que hay que importar.
- **`BREP.<guid>.smb`** — el **historial de diseño**: una instantánea del
  cuerpo tras cada feature del timeline (cientos de cuerpos intermedios, en
  coordenadas de modelado, sin las transformaciones aplicadas).

Los archivos antiguos (p.ej. la muestra del soporte de cortinas) sólo tienen
`.smb` y ahí las últimas instantáneas son los cuerpos finales. El importador
prefiere `.smbh` si existe. Los cuerpos ASM llevan además un puntero a un
record `transform` (3 ejes + traslación + escala) en `body[6]`; en las
muestras analizadas todos son identidad.

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

### 2.3 Offsets verificados de la topología

`values[0]` es siempre el `TypeName` de la entidad; los offsets siguientes
fueron verificados comprobando el **tipo** del record destino en todos los
records de cada clase:

| Record  | Posición → destino |
|---------|--------------------|
| body    | `[4]`→lump |
| lump    | `[5]`→shell, `[6]`→body |
| shell   | `[6]`→face, `[8]`→lump |
| face    | `[4]`→next face, `[5]`→loop, `[6]`→shell, `[8]`→surface |
| loop    | `[4]`→next loop (agujeros), `[5]`→coedge, `[6]`→face |
| coedge  | `[4]`→next, `[5]`→prev, `[7]`→edge, `[8]`=sense (enum), `[9]`→loop, `[11]`→pcurve |
| edge    | `[4]`→vértice inicio, `[5]`=t₀, `[6]`→vértice fin, `[7]`=t₁, `[9]`→curve, `[10]`=sentido (bool) |
| vertex  | `[6]`→point |

**Trampas descubiertas (importantes para cualquier lector):**

- La cadena `next` de coedges **no ordena** el recorrido del loop de forma
  fiable; hay que reconstruir el ciclo encadenando polilíneas por coincidencia
  de extremos (los vértices compartidos son exactos).
- El orden de los loops de una cara **no garantiza** el exterior primero; hay
  que ordenarlos por tamaño (bbox) — el exterior contiene a los agujeros.
- Una arista puede recorrer su curva **en contra**: entonces t₀/t₁ se
  refieren a la parametrización invertida `t → -t` y hay que **negarlos**
  antes de evaluar. Si no, los arcos circulares se muestrean en el
  **cuadrante especular** (arcos del rebaje "ojo de cerradura" del Perchero
  con 1.56 mm de error) y las intcurve extrapolan fuera de dominio (las
  paredes de las letras en relieve caían al fallback recto). El bool
  `edge[10]` **no es fiable** como flag de sentido: marca `False` algunas
  aristas invertidas, pero otras (intcurves con params fuera de dominio)
  llevan `True`. Regla robusta validada: muestrear con params directos y
  negados y quedarse con el candidato que **reproduce los vértices** de la
  arista (los vértices son exactos, el error del candidato correcto es 0).

### 2.4 Geometría de curvas

- **`straight`**: posiciones = punto base + dirección; params = distancia.
- **`ellipse`**: posiciones = centro, normal, eje mayor (|eje| = a);
  el primer escalar tras las 3 posiciones es el *ratio* menor/mayor.
  Los parámetros del edge (t₀, t₁) son **ángulos**:
  `P(t) = centro + eje·cos t + (n̂×êje)·(a·ratio)·sin t`.
- **`intcurve`/`exact_int_cur`, `spring_int_cur`, `par_int_cur`**: B-splines
  (ver §2.5). Las *spring* son los bordes de contacto de los fillets.

### 2.5 B-splines: `nubs` (no racional) y `nurbs` (racional)

Layout del bloque: `grado (int)`, `flag (int)`, `nº de nudos distintos (int)`,
luego pares `(valor double, multiplicidad int)`, y los **puntos de control**:
ternas (x,y,z) en `nubs`, cuádruplas (x,y,z,w) en `nurbs`.

Convenciones **no estándar** descubiertas (críticas):

- `flag=0` → curva abierta *clamped*, pero las multiplicidades de los nudos
  extremos se guardan como `grado` en vez de `grado+1`: hay que **sumar 1 a
  cada extremo** para reconstruir el vector de nudos. Validado: de Boor
  reproduce los extremos de todas las aristas spline con error 0.
- `flag=2` → curva **cerrada/periódica**. La reconstrucción *clamped* (+1 de
  multiplicidad en cada extremo) sigue siendo válida — la suma de
  multiplicidades coincide con los puntos de control almacenados y la curva
  cierra exacta (gap ~10⁻¹⁵) — pero las aristas sobre estas curvas pueden
  llevar **parámetros fuera del dominio** almacenado (p.ej. `[-T, 0]` con
  dominio `[0, T]`): hay que envolverlos **módulo el periodo** T =
  `nudo_final − nudo_inicial` antes de evaluar.
- El número de puntos de control **no se almacena**; se deriva de
  `len(nudos) − grado − 1` (evita leer de más hacia campos posteriores, como
  la tolerancia de ajuste).

### 2.6 Geometría de superficies

- **`plane`**: posiciones = origen, normal, referencia u.
- **`cone`** (incluye cilindros): posiciones = origen, eje, eje mayor
  (|eje mayor| = r₀); tras el ratio vienen `seno` y `coseno` del semiángulo.
  Parametrización: u = ángulo, v = distancia axial, `r(v) = r₀ + v·(sen/cos)`
  (el signo del taper se calibra empíricamente con el contorno).
  Validado: residuo < 6·10⁻⁶ en todos los puntos de contorno de la muestra.
- **`cyl_spl_sur`**: superficie de **extrusión** — una directriz B-spline
  barrida por un eje: `S(u,v) = C(u) + v·êje`. El eje es la tupla XYZ unitaria
  tras la definición; la directriz puede ser inline o `ref` (ver §2.7).
- **`torus`** (aparece en los `.smbh`, ASM 229): layout posicional
  `[6]=centro, [7]=eje, [8]=radio mayor R, [9]=radio menor r,
  [10]=dirección de referencia`. `S(u,v) = C + (R + r·cos v)·(cos u·M̂ +
  sin u·N̂) + r·sin v·êje`. R < r es válido (toro tipo limón/manzana en
  fillets de esquina). Validado: residuo < 1.4·10⁻⁴ en 102 caras.
  Ojo: el toro es periódico en **ambos** parámetros; una cara puede cruzar
  la costura de v (atan2 en ±π, el ecuador interior) igual que la de u, y el
  teselador debe desenvolver ambas (caras que la cruzaban salían 1.3–1.7 mm
  desplazadas al triangular un dominio uv roto).
- **`srf_srf_v_bl_spl_sur`**: **blend/fillet** entre dos superficies soporte
  (que van embebidas en el record como `blend_support_surface`, p.ej. un
  `cyl_spl_sur` y un `cone`). El record contiene la **curva del centro de la
  bola** como B-spline 3D a distancia constante = radio de ambos perfiles de
  contacto (curvas *spring*), más una B-spline **2D** periódica (ley del
  blend). Ojo: la superficie que Fusion evalúa/exporta **no** son arcos
  circulares de bola rodante — sus secciones transversales son mucho más
  planas (flecha ≈ 3 % de la cuerda en la muestra, frente al 21 % de un arco
  circular). La aproximación fiel al OBJ de referencia es la **franja reglada**
  entre puntos correspondientes de los dos perfiles de contacto.

### 2.7 Geometría compartida: el pool de `ref N`

ASM comparte datos de subtipo (curvas `*_int_cur`/`*_par_cur` y superficies
`*_spl_sur`). Regla de numeración (descubierta y validada geométricamente):

> Cada **definición** de subtipo — todo `TypeName` cuyo nombre termina en
> `_cur` o `_sur` — recibe el siguiente índice (base 0) en orden de aparición
> en el flujo, **excluyendo** los records `ATTRIB_CUSTOM`. Un valor `ref N`
> posterior referencia la definición N.

La primera vez que se usa un objeto compartido se serializa **inline**
(a menudo dentro de un `pcurve`); los usos siguientes emiten `ref N`.
Validación: 59 de 61 caras curvas de la muestra resuelven a superficies con
residuo < 10⁻³ (las 2 restantes son blends, sin evaluador cerrado).

### 2.8 `pcurve`

Curva en el espacio paramétrico (u,v) de una superficie (`exp_par_cur` +
`nubs` 2D). El primer `pcurve` sobre una superficie incluye la definición
completa de ésta (es su primera serialización; ver §2.7).

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
| Grafo topológico B-rep             | ✅ Offsets verificados (§2.3); cuerpos/caras/loops/aristas |
| Curvas (recta, elipse, nubs/nurbs) | ✅ de Boor validado; periódicas (flag=2) con wrap de parámetros |
| Superficies (plano, cono, toro, extrusión)| ✅ Evaluación analítica, residual < 6·10⁻⁶ (toro < 1.4·10⁻⁴) |
| Modelo vs. historial (`.smbh`/`.smb`) | ✅ Se importa el modelo actual; regla de ordinales con `delta_state`/`Begin` |
| Superficies blend (fillets)        | ✅ Franja reglada entre perfiles de contacto (§2.6) |
| Pool de referencias `ref N`        | ✅ Regla validada (§2.7) |
| Teselación de caras recortadas     | ✅ CDT en espacio paramétrico + loft (`tessellate.py`) |
| Fidelidad vs. OBJ de referencia    | ✅ Área 76.668 / 76.608 cm² (0.08 %); desviación media 10–26 µm, máx 140 µm; malla estanca (0 aristas no-manifold) |
| `BulkStream.dat` (árbol/materiales)| 🔜 Parcial (strings) |
| Operador de importación Blender    | ✅ Funcional (`importer.py`): escala, dedup de cuerpos, soldadura |
