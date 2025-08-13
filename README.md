# Sequal - Go ProForma Parser

A Go implementation of a ProForma 2.0 peptide sequence notation parser. This library provides comprehensive support for parsing, manipulating, and generating ProForma notation strings for peptide and protein sequences.

## Features

- **Complete ProForma 2.0 Support**: Parse and generate ProForma notation strings
- **Modification Handling**: Support for all modification types (static, ambiguous, labile, terminal)
- **Charge States**: Handle charge information and ionic species
- **Global Modifications**: Fixed modifications and isotope labeling
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
// Labile modification (typically glycans)
seq, _ := sequal.FromProforma("{Glycan:Hex(1)HexNAc(2)}PEPTIDE")
labileMods := seq.GetMods()[-3] // Labile modifications
fmt.Printf("Labile mod: %s\n", labileMods[0].GetValue()) // Output: Hex(1)HexNAc(2)
fmt.Printf("Mod type: %s\n", labileMods[0].GetModType()) // Output: labile
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
    "{InvalidGlycan}PEPTIDE",  // Invalid labile mod (must start with Glycan:)
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

## ProForma 2.0 Compliance

This implementation follows the ProForma 2.0 specification and supports:

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

## Contributing

Contributions are welcome! Please ensure all tests pass and add appropriate test coverage for new features.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

Based on the TypeScript [SequalJS](https://github.com/antmass/SequalJS) implementation. This Go port maintains compatibility while fixing some edge cases and providing Go-idiomatic APIs.