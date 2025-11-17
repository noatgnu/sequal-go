# Sequal - Go ProForma Parser

A Go implementation of a ProForma 2.1 peptide sequence notation parser. This library provides comprehensive support for parsing, manipulating, and generating ProForma notation strings for peptide and protein sequences.

## Features

- **Complete ProForma 2.1 Support**: Parse and generate ProForma notation strings with all 2.1 enhancements
- **ProForma 2.1 Features**:
  - Charged formulas with charge notation (e.g., `Formula:Zn1:z+2`)
  - Named entity definitions (e.g., `#g1:Glycan`)
  - Custom monosaccharide definitions
  - Terminal global modifications
  - Placement controls: Position constraints, limits, colocalization tags
  - Ion type notation for fragment ions
- **Modification Handling**: Support for all modification types (static, ambiguous, labile, terminal)
- **Charge States**: Handle charge information and ionic species
- **Global Modifications**: Fixed modifications and isotope labeling with placement controls
- **Complex Sequences**: Multi-chain, chimeric, and branched sequences
- **Error Handling**: Comprehensive validation and error reporting
- **Round-trip Conversion**: Parse ProForma → manipulate → regenerate ProForma

## Installation

```bash
go get github.com/noatgnu/sequal-go
```

## Quick Start

```go
package main

import (
    "fmt"
    "github.com/noatgnu/sequal-go/sequal"
)

func main() {
    // Parse a ProForma string
    seq, err := sequal.FromProforma("ELVIS[Phospho]K/2")
    if err != nil {
        panic(err)
    }

    // Access sequence information
    fmt.Printf("Sequence: %s\n", seq.ToStrippedString()) // Output: ELVISK
    fmt.Printf("Length: %d\n", seq.GetLength())          // Output: 6
    fmt.Printf("Charge: %d\n", *seq.GetCharge())         // Output: 2

    // Convert back to ProForma
    fmt.Printf("ProForma: %s\n", seq.ToProforma())       // Output: ELVIS[Phospho]K/2
}
```

## Basic Usage

### Parsing Simple Sequences

```go
// Basic peptide sequence
seq, _ := sequal.FromProforma("PEPTIDE")
fmt.Println(seq.ToStrippedString()) // Output: PEPTIDE

// Sequence with modification
seq, _ := sequal.FromProforma("PEP[+79.966]TIDE")
fmt.Println(seq.ToStrippedString()) // Output: PEPTIDE
```

### Working with Modifications

```go
// Named modification
seq, _ := sequal.FromProforma("ELVIS[Phospho]K")
aa := seq.GetSeq()[4] // Get the 'S' amino acid
mods := aa.GetMods()
fmt.Printf("Modification: %s\n", mods[0].GetValue()) // Output: Phospho

// Mass shift modification
seq, _ := sequal.FromProforma("PEP[+79.966]TIDE")
aa := seq.GetSeq()[2] // Get the 'P' amino acid
mods := aa.GetMods()
fmt.Printf("Mass shift: %f\n", *mods[0].GetMass()) // Output: 79.966000

// Multiple modifications on same residue
seq, _ := sequal.FromProforma("PEPS[Phospho][Acetyl]TIDE")
aa := seq.GetSeq()[3] // Get the 'S' amino acid
mods := aa.GetMods()
fmt.Printf("Mod 1: %s\n", mods[0].GetValue()) // Output: Phospho
fmt.Printf("Mod 2: %s\n", mods[1].GetValue()) // Output: Acetyl
```

### Terminal Modifications

```go
// N-terminal and C-terminal modifications
seq, _ := sequal.FromProforma("[Acetyl]-PEPTIDE-[Amidated]")

// Access terminal modifications
nTermMods := seq.GetMods()[-1] // N-terminal
cTermMods := seq.GetMods()[-2] // C-terminal

fmt.Printf("N-term: %s\n", nTermMods[0].GetValue()) // Output: Acetyl
fmt.Printf("C-term: %s\n", cTermMods[0].GetValue()) // Output: Amidated
```

### Ambiguous Modifications

```go
// Ambiguous modification (uncertain localization)
seq, _ := sequal.FromProforma("PEPS{Phospho}TIDE")
aa := seq.GetSeq()[3] // Get the 'S' amino acid
mods := aa.GetMods()
fmt.Printf("Mod type: %s\n", mods[0].GetModType()) // Output: ambiguous

// Ambiguous modification with localization score
seq, _ := sequal.FromProforma("PEPT[Phospho#1(0.75)]IDE")
aa := seq.GetSeq()[3] // Get the 'T' amino acid
mods := aa.GetMods()
fmt.Printf("Ambiguity group: %s\n", *mods[0].GetAmbiguityGroup()) // Output: 1
```

### Labile Modifications

```go
// Glycan labile modification
seq, _ := sequal.FromProforma("{Glycan:Hex(1)HexNAc(2)}PEPTIDE")
labileMods := seq.GetMods()[-3] // Labile modifications
fmt.Printf("Labile mod: %s\n", labileMods[0].GetValue()) // Output: Hex(1)HexNAc(2)
fmt.Printf("Mod type: %s\n", labileMods[0].GetModType()) // Output: labile

// Non-glycan labile modification
seq, _ := sequal.FromProforma("{Phospho}PEPTIDE")
labileMods = seq.GetMods()[-3]
fmt.Printf("Labile mod: %s\n", labileMods[0].GetValue()) // Output: Phospho

// Mass shift labile modification
seq, _ := sequal.FromProforma("{+79.966}PEPTIDE")
labileMods = seq.GetMods()[-3]
fmt.Printf("Labile mod: %s\n", labileMods[0].GetValue()) // Output: +79.966
```

### Charge States and Ionic Species

```go
// Simple charge state
seq, _ := sequal.FromProforma("PEPTIDE/2")
fmt.Printf("Charge: %d\n", *seq.GetCharge()) // Output: 2

// Negative charge
seq, _ := sequal.FromProforma("PEPTIDE/-3")
fmt.Printf("Charge: %d\n", *seq.GetCharge()) // Output: -3

// Charge with ionic species
seq, _ := sequal.FromProforma("PEPTIDE/2[+Na+]")
fmt.Printf("Charge: %d\n", *seq.GetCharge()) // Output: 2
fmt.Printf("Ionic species: %s\n", *seq.GetIonicSpecies()) // Output: +Na+
```

### Global Modifications

```go
// Fixed modification applied to all cysteines
seq, _ := sequal.FromProforma("<[Carbamidomethyl]@C>PEPTCDE")
globalMods := seq.GetGlobalMods()
fmt.Printf("Global mod: %s\n", globalMods[0].GetValue()) // Output: Carbamidomethyl

// Isotope labeling
seq, _ := sequal.FromProforma("<15N>PEPTIDE")
globalMods := seq.GetGlobalMods()
fmt.Printf("Isotope: %s\n", globalMods[0].GetValue()) // Output: 15N
```

### Unknown Position Modifications

```go
// Modification at unknown position
seq, _ := sequal.FromProforma("[Phospho]?PEPTIDE")
unknownMods := seq.GetMods()[-4] // Unknown position modifications
fmt.Printf("Unknown mod: %s\n", unknownMods[0].GetValue()) // Output: Phospho

// Multiple unknown modifications
seq, _ := sequal.FromProforma("[Phospho]^2?PEPTIDE")
unknownMods := seq.GetMods()[-4]
fmt.Printf("Count: %d\n", len(unknownMods)) // Output: 2
```

### Sequence Ambiguity

```go
// Ambiguous amino acid sequence
seq, _ := sequal.FromProforma("(?LI)PEPTIDESEQUENCE")
ambiguities := seq.GetSequenceAmbiguities()
fmt.Printf("Ambiguous sequence: %s\n", ambiguities[0].GetValue()) // Output: LI
fmt.Printf("Position: %d\n", ambiguities[0].GetPosition()) // Output: 0 (always 0)
```

### Gap Notation

```go
// Gap with mass information
seq, _ := sequal.FromProforma("RTAAX[+367.0537]WT")
aa := seq.GetSeq()[4] // Get the 'X' amino acid
mods := aa.GetMods()
fmt.Printf("Gap mass: %f\n", *mods[0].GetMass()) // Output: 367.053700
fmt.Printf("Mod type: %s\n", mods[0].GetModType()) // Output: gap
```

### Multi-Chain Sequences

```go
// Multiple peptide chains
seq, _ := sequal.FromProforma("PEPTIDE//SEQUENCE//THIRD")
fmt.Printf("Is multi-chain: %t\n", seq.IsMultiChain()) // Output: true

chains := seq.GetChains()
fmt.Printf("Chain count: %d\n", len(chains)) // Output: 3

for i, chain := range chains {
    fmt.Printf("Chain %d: %s\n", i+1, chain.ToStrippedString())
}
// Output:
// Chain 1: PEPTIDE
// Chain 2: SEQUENCE  
// Chain 3: THIRD
```

### Chimeric Sequences

```go
// Chimeric peptide (multiple peptidoforms)
seq, _ := sequal.FromProforma("PEPTIDE/2+ANOTHER/3")
fmt.Printf("Is chimeric: %t\n", seq.IsChimeric()) // Output: true

peptidoforms := seq.GetPeptidoforms()
fmt.Printf("Peptidoform count: %d\n", len(peptidoforms)) // Output: 2

for i, pep := range peptidoforms {
    fmt.Printf("Peptidoform %d: %s\n", i+1, pep.ToStrippedString())
    if charge := pep.GetCharge(); charge != nil {
        fmt.Printf("  Charge: %d\n", *charge)
    }
}
```

## ProForma 2.1 Features

### Charged Formulas

```go
// Charged formula with charge notation
seq, _ := sequal.FromProforma("PEPT[Formula:Zn1:z+2]IDE")
aa := seq.GetSeq()[3] // Get the 'T' amino acid
mods := aa.GetMods()
modValue := mods[0].GetModificationValue()

// Access charged formula information
for _, pv := range modValue.GetPipeValues() {
    if pv.GetCharge() != nil {
        fmt.Printf("Charge: %s\n", *pv.GetCharge()) // Output: z+2
        fmt.Printf("Charge value: %d\n", *pv.GetChargeValue()) // Output: 2
    }
}

// Multiple charged formulas
seq, _ := sequal.FromProforma("PE[Formula:Cu1:z+1]PT[Formula:Zn1:z+2]IDE")
```

### Placement Controls

```go
// Position constraint on global modification
seq, _ := sequal.FromProforma("<[Oxidation|Position:M]@M>PEPTMIDE")
globalMod := seq.GetGlobalMods()[0]
positions := globalMod.GetPositionConstraint()
fmt.Printf("Positions: %v\n", positions) // Output: [M]

// Multiple position constraints
seq, _ := sequal.FromProforma("<[Phospho|Position:S,T,Y]@S,T,Y>PEPTIDES")
globalMod = seq.GetGlobalMods()[0]
positions = globalMod.GetPositionConstraint()
fmt.Printf("Positions: %v\n", positions) // Output: [S T Y]

// Limit per position
seq, _ := sequal.FromProforma("<[Phospho|Limit:2]@S,T,Y>STYSTY")
globalMod = seq.GetGlobalMods()[0]
limit := globalMod.GetLimitPerPosition()
fmt.Printf("Limit: %d\n", *limit) // Output: 2

// Colocalize modifications of known position (CoMKP)
seq, _ := sequal.FromProforma("<[Oxidation|CoMKP]@M>PEPTIDE")
globalMod = seq.GetGlobalMods()[0]
fmt.Printf("CoMKP: %t\n", globalMod.GetColocalizeKnown()) // Output: true

// Colocalize modifications of unknown position (CoMUP)
seq, _ := sequal.FromProforma("<[Oxidation|CoMUP]@M>PEPTIDE")
globalMod = seq.GetGlobalMods()[0]
fmt.Printf("CoMUP: %t\n", globalMod.GetColocalizeUnknown()) // Output: true

// Combined placement controls
seq, _ := sequal.FromProforma("<[Phospho|Position:S,T,Y|Limit:2|CoMKP]@S,T,Y>PEPTIDES")
globalMod = seq.GetGlobalMods()[0]
fmt.Printf("Positions: %v\n", globalMod.GetPositionConstraint()) // Output: [S T Y]
fmt.Printf("Limit: %d\n", *globalMod.GetLimitPerPosition()) // Output: 2
fmt.Printf("CoMKP: %t\n", globalMod.GetColocalizeKnown()) // Output: true
```

### Ion Type Notation

```go
// Fragment ion types
seq, _ := sequal.FromProforma("PEPT[a-type-ion]IDE")
aa := seq.GetSeq()[3] // Get the 'T' amino acid
mods := aa.GetMods()
fmt.Printf("Is ion type: %t\n", mods[0].IsIonType()) // Output: true

// All supported ion types: a, b, c, x, y, z
ionTypes := []string{
    "PEPT[a-type-ion]IDE",
    "PEPT[b-type-ion]IDE",
    "PEPT[c-type-ion]IDE",
    "PEPT[x-type-ion]IDE",
    "PEPT[y-type-ion]IDE",
    "PEPT[z-type-ion]IDE",
}

// Using Unimod ion references
seq, _ := sequal.FromProforma("PEPT[UNIMOD:140]IDE") // a-type-ion
aa = seq.GetSeq()[3]
mods = aa.GetMods()
fmt.Printf("Is ion type: %t\n", mods[0].IsIonType()) // Output: true

// Short form Unimod reference
seq, _ := sequal.FromProforma("PEPT[U:2132]IDE") // b-type-ion
aa = seq.GetSeq()[3]
mods = aa.GetMods()
fmt.Printf("Is ion type: %t\n", mods[0].IsIonType()) // Output: true

// Case insensitive
seq, _ := sequal.FromProforma("PEPT[A-TYPE-ION]IDE")
aa = seq.GetSeq()[3]
mods = aa.GetMods()
fmt.Printf("Is ion type: %t\n", mods[0].IsIonType()) // Output: true
```

### Complex ProForma 2.1 Combinations

```go
// Placement controls with ion notation
seq, _ := sequal.FromProforma("<[Phospho|Position:S,T,Y|Limit:2]@S,T,Y>PEPT[a-type-ion]IDE")

// Charged formulas with ion types
seq, _ := sequal.FromProforma("PE[a-type-ion]PT[Formula:Zn1:z+2]IDE")

// Multiple global modifications with placement controls
seq, _ := sequal.FromProforma("<[Phospho|Position:S,T,Y|CoMKP]@S,T,Y><[Oxidation|Position:M|CoMUP]@M>STMYST")

// Full feature combination
seq, _ := sequal.FromProforma("<[Phospho|Position:S,T,Y|Limit:2|CoMKP]@S,T,Y>[Acetyl]-PE[a-type-ion]PT[Formula:Zn1:z+2]IDE")
```

## Advanced Features

### Crosslinks and Branches

```go
// Crosslinked peptides
seq, _ := sequal.FromProforma("PEPTK[XL:DSS#XL1]IDE")
aa := seq.GetSeq()[4] // Get the 'K' amino acid
mods := aa.GetMods()
fmt.Printf("Crosslink ID: %s\n", *mods[0].GetCrosslinkID()) // Output: XL1

// Branch modifications
seq, _ := sequal.FromProforma("PEPTK[Branch#BRANCH]IDE")
aa := seq.GetSeq()[4] // Get the 'K' amino acid
mods := aa.GetMods()
fmt.Printf("Mod type: %s\n", mods[0].GetModType()) // Output: branch
```

### Information Tags

```go
// Modifications with information tags
seq, _ := sequal.FromProforma("ELVIS[Phospho|INFO:newly discovered]K")
aa := seq.GetSeq()[4] // Get the 'S' amino acid
mods := aa.GetMods()
infoTags := mods[0].GetInfoTags()
fmt.Printf("Info tag: %s\n", infoTags[0]) // Output: newly discovered
```

### Range Modifications

```go
// Modification applied to a range of residues
seq, _ := sequal.FromProforma("(PEP)[+79.966]TIDE")
fmt.Printf("Sequence: %s\n", seq.ToStrippedString()) // Output: PEPTIDE

// Check modifications on range positions (0, 1, 2 for P, E, P)
for i := 0; i < 3; i++ {
    aa := seq.GetSeq()[i]
    if len(aa.GetMods()) > 0 {
        fmt.Printf("Position %d has modification: %s\n", i, aa.GetMods()[0].GetValue())
    }
}
```

## Utility Functions

### Sequence Analysis

```go
seq, _ := sequal.FromProforma("PEPTIDE")

// Count specific amino acids
count := seq.Count("P", 0, seq.GetLength())
fmt.Printf("P count: %d\n", count) // Output: 2

// Find gaps in sequence
gapSeq, _ := sequal.FromProforma("PEP-TIDE")
gaps := gapSeq.Gaps()
fmt.Printf("Has gaps: %v\n", gaps) // Output: [false false false true false false false false]

// Convert to map representation
seqMap := seq.ToMap()
fmt.Printf("Map: %+v\n", seqMap)
```

### Pattern Matching

```go
seq, _ := sequal.FromProforma("PEPTIDES")

// Find amino acid positions using regex
positions, _ := seq.FindWithRegex("S", nil)
if len(positions) > 0 {
    fmt.Printf("'S' found at position: %d\n", positions[0][0]) // Output: 7
}

// Find patterns
positions, _ = seq.FindWithRegex("P.P", nil) // Find P followed by any char then P
fmt.Printf("Pattern matches: %v\n", positions)
```

## Error Handling

The library provides comprehensive error handling for malformed ProForma strings:

```go
// Various error cases
testCases := []string{
    "ELVIS[PhosphoPEPTIDE",    // Unclosed bracket
    "<15NPEPTIDE",             // Unclosed global mod
    "ELVIS(PEPTIDE",           // Unclosed parenthesis
    "ELVIS)PEPTIDE",           // Unmatched closing parenthesis
}

for _, test := range testCases {
    if _, err := sequal.FromProforma(test); err != nil {
        fmt.Printf("Error for '%s': %s\n", test, err)
    }
}
```

## Round-trip Conversion

The library supports parsing ProForma strings and regenerating them:

```go
original := "ELVIS[Phospho]K/2[+Na+]"
seq, _ := sequal.FromProforma(original)
regenerated := seq.ToProforma()
fmt.Printf("Original:    %s\n", original)
fmt.Printf("Regenerated: %s\n", regenerated)
// Both should be identical for valid ProForma strings
```

## Testing

Run the comprehensive test suite:

```bash
go test ./sequal -v
```

The tests cover:
- Basic sequence parsing
- All modification types
- Charge states and ionic species
- Global modifications
- Multi-chain and chimeric sequences
- Error handling
- Round-trip conversion
- Complex ProForma examples
- ProForma 2.1 features:
  - Charged formulas
  - Placement controls (Position, Limit, CoMKP, CoMUP)
  - Ion type notation
  - Integration tests combining multiple 2.1 features

## ProForma 2.1 Compliance

This implementation follows the ProForma 2.1 specification and supports all ProForma 2.0 features plus new 2.1 enhancements:

### ProForma 2.0 Features
- Basic amino acid sequences
- Static modifications `[mod]`
- Ambiguous modifications `{mod}`
- Terminal modifications `[mod]-` and `-[mod]`
- Labile modifications `{Glycan:mod}`
- Global modifications `<mod>`
- Mass shifts `[+/-mass]`
- Charge states `/charge`
- Ionic species `/charge[species]`
- Unknown position mods `[mod]?`
- Sequence ambiguity `(?seq)`
- Range modifications `(range)[mod]`
- Crosslinks `[XL:name#ID]`
- Branches `[Branch#ID]`
- Information tags `[mod|INFO:info]`
- Multi-chain sequences `seq//seq`
- Chimeric sequences `seq+seq`

### ProForma 2.1 Enhancements
- **Charged Formulas** (Section 11.1): `[Formula:Zn1:z+2]` - Formulas with charge notation
- **Named Entity Definitions** (Section 11.3): `[#g1:Glycan]` - Named entity references
- **Custom Monosaccharide Definitions** (Section 11.4): Custom glycan building blocks
- **Terminal Global Modifications** (Section 11.5): Global modifications on termini
- **Placement Controls** (Section 11.2): Advanced global modification constraints
  - Position constraints: `[mod|Position:M]@M`
  - Limit per position: `[mod|Limit:2]@S,T,Y`
  - CoMKP (Colocalize with known positions): `[mod|CoMKP]@M`
  - CoMUP (Colocalize with unknown positions): `[mod|CoMUP]@M`
- **Ion Notation** (Section 11.6): Fragment ion type semantic flag
  - Named ion types: `[a-type-ion]`, `[b-type-ion]`, `[c-type-ion]`, `[x-type-ion]`, `[y-type-ion]`, `[z-type-ion]`
  - Unimod ion references: `[UNIMOD:140]`, `[U:2132]`

## Contributing

Contributions are welcome! Please ensure all tests pass and add appropriate test coverage for new features.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

Based on the TypeScript [SequalJS](https://github.com/antmass/SequalJS) implementation. This Go port maintains compatibility while fixing some edge cases and providing Go-idiomatic APIs.