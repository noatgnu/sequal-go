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
	"U": 255.15829, // Note: U appears twice in original Python code
}

// GlycanBlockDict maps glycan block names to their masses
var GlycanBlockDict = map[string]float64{
	"HexNAc":  203.079372520,
	"Hex":     162.0528234185,
	"Fuc":     146.057908799,
	"NeuAc":   291.0954165066,
	"Sulfo":   79.9568148602,
	"Phospho": 79.9663305228,
	"Pent":    132.0422587348,
	"NeuGc":   307.0903311261,
}

// Monosaccharides is a set of known monosaccharide names
// In Go, we represent sets as maps with bool values
var Monosaccharides = map[string]bool{
	"Hex":     true,
	"HexNAc":  true,
	"HexS":    true,
	"HexP":    true,
	"HexNAcS": true,
	"dHex":    true,
	"NeuAc":   true,
	"NeuGc":   true,
	"Pen":     true,
	"Fuc":     true,
}
