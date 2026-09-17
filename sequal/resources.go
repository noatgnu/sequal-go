// Package sequal provides a library for parsing and generating ProForma 2.1 strings.
// This file contains constants and data maps used throughout the library.
package sequal

// Atomic mass constants
const (
	Proton = 1.007277
	H      = 1.007825
	O      = 15.99491463
)

// AAMass maps amino acid one-letter codes to their masses
var AAMass = map[string]float64{
	"A": 71.037114,
	"R": 156.101111,
	"N": 114.042927,
	"D": 115.026943,
	"C": 103.009185,
	"E": 129.042593,
	"Q": 128.058578,
	"G": 57.021464,
	"H": 137.058912,
	"I": 113.084064,
	"L": 113.084064,
	"K": 128.094963,
	"M": 131.040485,
	"F": 147.068414,
	"P": 97.052764,
	"S": 87.032028,
	"T": 101.047679,
	"W": 186.079313,
	"Y": 163.06332,
	"V": 99.068414,
	"X": 0,
	"O": 150.03794,
	"U": 255.15829,  // Note: U appears twice in original Python code
	"B": 114.534935, // avg of Asp/Asn
	"Z": 128.550586, // avg of Glu/Gln
	"J": 113.084064, // Leu/Ile, exact
}

// GlycanBlockDict maps glycan block names to their monoisotopic residue masses (Da)
var GlycanBlockDict = map[string]float64{
	"Hex":     162.0528234185,
	"HexNAc":  203.079372520,
	"HexS":    242.009638,
	"HexP":    242.019154,
	"HexNAcS": 283.036187,
	"HexN":    161.068808,
	"HexNS":   241.025623,
	"dHex":    146.057908799,
	"Fuc":     146.057908799,
	"aHex":    176.032088,
	"en,aHex": 158.021523,
	"Neu":     249.084852,
	"NeuAc":   291.0954165066,
	"NeuGc":   307.0903311261,
	"Sug":     42.010565,
	"Tri":     72.021129,
	"Tet":     102.031694,
	"Pen":     132.0422587348,
	"Hep":     192.063388,
	"Oct":     222.073953,
	"Non":     252.084517,
	"Dec":     282.095082,
	"Kdo":     220.058303,
	"Kdn":     250.068867,
	"Sulfo":   79.9568148602,
	"Phospho": 79.9663305228,
}

// Monosaccharides is a set of known monosaccharide names
// In Go, we represent sets as maps with bool values
var Monosaccharides = map[string]bool{
	"Hex":     true,
	"HexNAc":  true,
	"HexS":    true,
	"HexP":    true,
	"HexNAcS": true,
	"HexN":    true,
	"HexNS":   true,
	"dHex":    true,
	"Fuc":     true,
	"aHex":    true,
	"en,aHex": true,
	"Neu":     true,
	"NeuAc":   true,
	"NeuGc":   true,
	"Sug":     true,
	"Tri":     true,
	"Tet":     true,
	"Pen":     true,
	"Hep":     true,
	"Oct":     true,
	"Non":     true,
	"Dec":     true,
	"Kdo":     true,
	"Kdn":     true,
	"Sulfo":   true,
	"Phospho": true,
}
