# Contrato de datos

El pipeline escribe y la web lee `web/data/<familia>.json`.

```ts
interface FamilyData {
  family: string                 // "GH1"
  source: "pipeline"
  generated: string              // fecha ISO
  reference_pdb: string
  ligand: string
  reference_positions: number[]  // [165, 360]
  fits: Record<"catalytic" | "core", {
    // Columnas paralelas: una entrada por residuo
    pdb: string[]                // código completo, 4 letras
    chain: string[]              // "" si la fuente no lo guarda
    auth_seq_id: (number | null)[]
    hmm_position: (number | null)[]   // null si no tiene posición
    residue: string[]            // ALA, TRP...
    domain: string[]             // Eukaryota, Bacteria, Archaea, Unknown
    theta: number[]              // radianes, desde +y
    phi: number[]                // radianes, atan2(z, x), de 0 a 2π
    distance: number[]           // Å al origen
    x: number[]; y: number[]; z: number[]
  } | null>
  positions: { hmm_position: number; n: number; rmsf: number; distance: number; top_residue: string }[]
  horizon: { grid: number; north: number[]; south: number[] } | null   // mediana en Å por celda
}
```

Sistema de referencia: origen en el carbono anomérico del azúcar −1; eje y hacia el oxígeno glicosídico; eje x hacia la proyección del C4 en y = 0.

Ficheros: todos llevan `source: "pipeline"` (D-023).

Vista: desde el centro. Norte antihorario, sur horario, φ = 0 arriba.
