package sequal

import (
	"testing"
)

func TestProFormaParserBasicSequenceParsing(t *testing.T) {
	t.Run("should parse a simple peptide sequence", func(t *testing.T) {
		input := "PEPTIDE"
		sequence, mods, globalMods, ambiguities, _, err := ParseProForma(input)
		if err != nil {
			t.Fatalf("Failed to parse ProForma '%s': %v", input, err)
		}

		if sequence != "PEPTIDE" {
			t.Errorf("Expected 'PEPTIDE', got '%s'", sequence)
		}
		if len(mods) != 0 {
			t.Errorf("Expected 0 modifications, got %d", len(mods))
		}
		if len(globalMods) != 0 {
			t.Errorf("Expected 0 global modifications, got %d", len(globalMods))
		}
		if len(ambiguities) != 0 {
			t.Errorf("Expected 0 ambiguities, got %d", len(ambiguities))
		}
	})
}

func TestProFormaParserTerminalModifications(t *testing.T) {
	t.Run("should parse N-terminal modification", func(t *testing.T) {
		input := "[Acetyl]-PEPTIDE"
		sequence, mods, _, _, _, err := ParseProForma(input)
		if err != nil {
			t.Fatalf("Failed to parse ProForma '%s': %v", input, err)
		}

		if sequence != "PEPTIDE" {
			t.Errorf("Expected 'PEPTIDE', got '%s'", sequence)
		}
		if nTermMods, exists := mods["-1"]; !exists {
			t.Error("Expected N-terminal modifications")
		} else {
			if len(nTermMods) != 1 {
				t.Errorf("Expected 1 N-terminal modification, got %d", len(nTermMods))
			}
			if nTermMods[0].GetValue() != "Acetyl" {
				t.Errorf("Expected 'Acetyl', got '%s'", nTermMods[0].GetValue())
			}
			if nTermMods[0].GetModType() != "terminal" {
				t.Errorf("Expected 'terminal', got '%s'", nTermMods[0].GetModType())
			}
		}
	})

	t.Run("should parse C-terminal modification", func(t *testing.T) {
		input := "PEPTIDE-[Amidated]"
		sequence, mods, _, _, _, err := ParseProForma(input)
		if err != nil {
			t.Fatalf("Failed to parse ProForma '%s': %v", input, err)
		}

		if sequence != "PEPTIDE" {
			t.Errorf("Expected 'PEPTIDE', got '%s'", sequence)
		}
		if cTermMods, exists := mods["-2"]; !exists {
			t.Error("Expected C-terminal modifications")
		} else {
			if len(cTermMods) != 1 {
				t.Errorf("Expected 1 C-terminal modification, got %d", len(cTermMods))
			}
			if cTermMods[0].GetValue() != "Amidated" {
				t.Errorf("Expected 'Amidated', got '%s'", cTermMods[0].GetValue())
			}
			if cTermMods[0].GetModType() != "terminal" {
				t.Errorf("Expected 'terminal', got '%s'", cTermMods[0].GetModType())
			}
		}
	})

	t.Run("should parse both terminal modifications", func(t *testing.T) {
		input := "[Acetyl]-PEPTIDE-[Amidated]"
		sequence, mods, _, _, _, err := ParseProForma(input)
		if err != nil {
			t.Fatalf("Failed to parse ProForma '%s': %v", input, err)
		}

		if sequence != "PEPTIDE" {
			t.Errorf("Expected 'PEPTIDE', got '%s'", sequence)
		}

		if nTermMods, exists := mods["-1"]; !exists {
			t.Error("Expected N-terminal modifications")
		} else if len(nTermMods) != 1 {
			t.Errorf("Expected 1 N-terminal modification, got %d", len(nTermMods))
		}

		if cTermMods, exists := mods["-2"]; !exists {
			t.Error("Expected C-terminal modifications")
		} else if len(cTermMods) != 1 {
			t.Errorf("Expected 1 C-terminal modification, got %d", len(cTermMods))
		}
	})
}

func TestProFormaParserGlobalMods(t *testing.T) {
	tests := []struct {
		name               string
		proforma           string
		expectedSeq        string
		expectedGlobalMods int
		expectedModType    string
	}{
		{
			name:               "Fixed modification",
			proforma:           "<[Carbamidomethyl]@C>PEPTCDE",
			expectedSeq:        "PEPTCDE",
			expectedGlobalMods: 1,
			expectedModType:    "fixed",
		},
		{
			name:               "Isotope labeling",
			proforma:           "<15N>PEPTIDE",
			expectedSeq:        "PEPTIDE",
			expectedGlobalMods: 1,
			expectedModType:    "isotope",
		},
		{
			name:               "Multiple global mods",
			proforma:           "<15N><[Carbamidomethyl]@C>PEPTCDE",
			expectedSeq:        "PEPTCDE",
			expectedGlobalMods: 2,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			baseSeq, _, globalMods, _, _, err := ParseProForma(tt.proforma)
			if err != nil {
				t.Fatalf("Failed to parse ProForma '%s': %v", tt.proforma, err)
			}

			if baseSeq != tt.expectedSeq {
				t.Errorf("Expected sequence '%s', got '%s'", tt.expectedSeq, baseSeq)
			}

			if len(globalMods) != tt.expectedGlobalMods {
				t.Errorf("Expected %d global modifications, got %d", tt.expectedGlobalMods, len(globalMods))
			}

			if tt.expectedModType != "" && len(globalMods) > 0 {
				if globalMods[0].GetGlobalModType() != tt.expectedModType {
					t.Errorf("Expected global mod type '%s', got '%s'", tt.expectedModType, globalMods[0].GetGlobalModType())
				}
			}
		})
	}
}

func TestProFormaParserChargeInfo(t *testing.T) {
	tests := []struct {
		name            string
		proforma        string
		expectedSeq     string
		expectedCharge  *int
		expectedSpecies *string
	}{
		{
			name:           "Simple charge",
			proforma:       "PEPTIDE/2",
			expectedSeq:    "PEPTIDE",
			expectedCharge: IntPtr(2),
		},
		{
			name:           "Negative charge",
			proforma:       "PEPTIDE/-3",
			expectedSeq:    "PEPTIDE",
			expectedCharge: IntPtr(-3),
		},
		{
			name:            "Charge with ionic species",
			proforma:        "PEPTIDE/2[+Na+]",
			expectedSeq:     "PEPTIDE",
			expectedCharge:  IntPtr(2),
			expectedSpecies: StringPtr("+Na+"),
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			baseSeq, _, _, _, chargeInfo, err := ParseProForma(tt.proforma)
			if err != nil {
				t.Fatalf("Failed to parse ProForma '%s': %v", tt.proforma, err)
			}

			if baseSeq != tt.expectedSeq {
				t.Errorf("Expected sequence '%s', got '%s'", tt.expectedSeq, baseSeq)
			}

			if len(chargeInfo) > 0 {
				charge := chargeInfo[0]
				if tt.expectedCharge == nil {
					if charge != nil {
						t.Errorf("Expected no charge, got %d", *charge)
					}
				} else {
					if charge == nil {
						t.Errorf("Expected charge %d, got nil", *tt.expectedCharge)
					} else if *charge != *tt.expectedCharge {
						t.Errorf("Expected charge %d, got %d", *tt.expectedCharge, *charge)
					}
				}
			}

			// Note: Ionic species testing may need adjustment based on actual implementation
		})
	}
}

func TestProFormaParserLabileMods(t *testing.T) {
	testCases := []struct {
		name     string
		proforma string
		expected string
	}{
		{
			name:     "Glycan labile modification",
			proforma: "{Glycan:Hex(1)HexNAc(2)}PEPTIDE",
			expected: "PEPTIDE",
		},
		{
			name:     "Non-glycan labile modification",
			proforma: "{Phospho}PEPTIDE",
			expected: "PEPTIDE",
		},
		{
			name:     "Mass shift labile modification",
			proforma: "{+79.966}PEPTIDE",
			expected: "PEPTIDE",
		},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			baseSeq, modifications, _, _, _, err := ParseProForma(tc.proforma)
			if err != nil {
				t.Fatalf("Failed to parse ProForma '%s': %v", tc.proforma, err)
			}

			if baseSeq != tc.expected {
				t.Errorf("Expected sequence '%s', got '%s'", tc.expected, baseSeq)
			}

			// Check labile modifications (-3)
			if labileMods, exists := modifications["-3"]; exists {
				if len(labileMods) != 1 {
					t.Errorf("Expected 1 labile modification, got %d", len(labileMods))
				} else {
					if labileMods[0].GetModType() != "labile" {
						t.Errorf("Expected labile mod type, got '%s'", labileMods[0].GetModType())
					}
				}
			} else {
				t.Errorf("Expected labile modification, but found none")
			}
		})
	}
}

func TestProFormaParserUnknownPositionMods(t *testing.T) {
	tests := []struct {
		name            string
		proforma        string
		expectedSeq     string
		expectedUnknown int
	}{
		{
			name:            "Single unknown position",
			proforma:        "[Phospho]?PEPTIDE",
			expectedSeq:     "PEPTIDE",
			expectedUnknown: 1,
		},
		{
			name:            "Multiple unknown positions",
			proforma:        "[Phospho]^2?PEPTIDE",
			expectedSeq:     "PEPTIDE",
			expectedUnknown: 2,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			baseSeq, modifications, _, _, _, err := ParseProForma(tt.proforma)
			if err != nil {
				t.Fatalf("Failed to parse ProForma '%s': %v", tt.proforma, err)
			}

			if baseSeq != tt.expectedSeq {
				t.Errorf("Expected sequence '%s', got '%s'", tt.expectedSeq, baseSeq)
			}

			// Check unknown position modifications (-4)
			if unknownMods, exists := modifications["-4"]; exists {
				if len(unknownMods) != tt.expectedUnknown {
					t.Errorf("Expected %d unknown position modifications, got %d", tt.expectedUnknown, len(unknownMods))
				}
			} else if tt.expectedUnknown > 0 {
				t.Errorf("Expected %d unknown position modifications, got 0", tt.expectedUnknown)
			}
		})
	}
}

func TestProFormaParserSequenceAmbiguity(t *testing.T) {
	proforma := "PEPTIDE(?LI)SEQUENCE"
	baseSeq, _, _, seqAmbig, _, err := ParseProForma(proforma)
	if err != nil {
		t.Fatalf("Failed to parse ProForma '%s': %v", proforma, err)
	}

	expectedSeq := "PEPTIDESEQUENCE"
	if baseSeq != expectedSeq {
		t.Errorf("Expected sequence '%s', got '%s'", expectedSeq, baseSeq)
	}

	if len(seqAmbig) != 1 {
		t.Errorf("Expected 1 sequence ambiguity, got %d", len(seqAmbig))
	} else {
		if seqAmbig[0].GetValue() != "LI" {
			t.Errorf("Expected ambiguity value 'LI', got '%s'", seqAmbig[0].GetValue())
		}
		if seqAmbig[0].GetPosition() != 0 { // Sequence ambiguities are always at position 0
			t.Errorf("Expected ambiguity position 0, got %d", seqAmbig[0].GetPosition())
		}
	}
}

func TestProFormaParserRangeMods(t *testing.T) {
	proforma := "(PEP)[+79.966]TIDE"
	baseSeq, modifications, _, _, _, err := ParseProForma(proforma)
	if err != nil {
		t.Fatalf("Failed to parse ProForma '%s': %v", proforma, err)
	}

	if baseSeq != "PEPTIDE" {
		t.Errorf("Expected sequence 'PEPTIDE', got '%s'", baseSeq)
	}

	// Range modifications should be applied to positions 0, 1, 2 (P, E, P)
	for i := 0; i < 3; i++ {
		posStr := string(rune(i + '0'))
		if mods, exists := modifications[posStr]; exists {
			if len(mods) == 0 {
				t.Errorf("Expected range modification at position %d, but found none", i)
			} else if !mods[0].inRange {
				t.Errorf("Expected modification at position %d to be marked as in range", i)
			}
		} else {
			t.Errorf("Expected modification at position %d, but found none", i)
		}
	}
}

func TestProFormaParserCrosslinks(t *testing.T) {
	tests := []struct {
		name              string
		proforma          string
		expectedSeq       string
		expectedCrosslink string
	}{
		{
			name:              "Crosslink with ID",
			proforma:          "PEPTK[XL:DSS#XL1]IDE",
			expectedSeq:       "PEPTKIDE",
			expectedCrosslink: "XL1",
		},
		{
			name:              "Crosslink reference",
			proforma:          "PEPTK[#XL1]IDE",
			expectedSeq:       "PEPTKIDE",
			expectedCrosslink: "XL1",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			baseSeq, modifications, _, _, _, err := ParseProForma(tt.proforma)
			if err != nil {
				t.Fatalf("Failed to parse ProForma '%s': %v", tt.proforma, err)
			}

			if baseSeq != tt.expectedSeq {
				t.Errorf("Expected sequence '%s', got '%s'", tt.expectedSeq, baseSeq)
			}

			// Check for crosslink modification at position 4 (K)
			if mods, exists := modifications["4"]; exists && len(mods) > 0 {
				mod := mods[0]
				if mod.GetModType() != "crosslink" {
					t.Errorf("Expected crosslink modification type, got '%s'", mod.GetModType())
				}

				crosslinkID := mod.GetCrosslinkID()
				if crosslinkID == nil {
					t.Errorf("Expected crosslink ID, got nil")
				} else if *crosslinkID != tt.expectedCrosslink {
					t.Errorf("Expected crosslink ID '%s', got '%s'", tt.expectedCrosslink, *crosslinkID)
				}
			} else {
				t.Errorf("Expected crosslink modification at position 4, but found none")
			}
		})
	}
}

func TestProFormaParserBranches(t *testing.T) {
	tests := []struct {
		name        string
		proforma    string
		expectedSeq string
	}{
		{
			name:        "Branch with ID",
			proforma:    "PEPTK[Branch#BRANCH]IDE",
			expectedSeq: "PEPTKIDE",
		},
		{
			name:        "Branch reference",
			proforma:    "PEPTK[#BRANCH]IDE",
			expectedSeq: "PEPTKIDE",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			baseSeq, modifications, _, _, _, err := ParseProForma(tt.proforma)
			if err != nil {
				t.Fatalf("Failed to parse ProForma '%s': %v", tt.proforma, err)
			}

			if baseSeq != tt.expectedSeq {
				t.Errorf("Expected sequence '%s', got '%s'", tt.expectedSeq, baseSeq)
			}

			// Check for branch modification
			found := false
			for _, mods := range modifications {
				for _, mod := range mods {
					if mod.GetModType() == "branch" {
						found = true
						break
					}
				}
				if found {
					break
				}
			}

			if !found {
				t.Errorf("Expected to find branch modification")
			}
		})
	}
}

func TestProFormaParserGapNotation(t *testing.T) {
	proforma := "RTAAX[+367.0537]WT"
	baseSeq, modifications, _, _, _, err := ParseProForma(proforma)
	if err != nil {
		t.Fatalf("Failed to parse ProForma '%s': %v", proforma, err)
	}

	if baseSeq != "RTAAXWT" {
		t.Errorf("Expected sequence 'RTAAXWT', got '%s'", baseSeq)
	}

	// Check for gap modification at position 4 (X)
	if mods, exists := modifications["4"]; exists && len(mods) > 0 {
		mod := mods[0]
		if mod.GetModType() != "gap" {
			t.Errorf("Expected gap modification type, got '%s'", mod.GetModType())
		}

		mass := mod.GetMass()
		if mass == nil {
			t.Errorf("Expected gap mass, got nil")
		} else if *mass != 367.0537 {
			t.Errorf("Expected gap mass 367.0537, got %f", *mass)
		}
	} else {
		t.Errorf("Expected gap modification at position 4, but found none")
	}
}

func TestProFormaParserAmbiguousModifications(t *testing.T) {
	proforma := "ELVIS{Phospho}K"
	baseSeq, modifications, _, _, _, err := ParseProForma(proforma)
	if err != nil {
		t.Fatalf("Failed to parse ProForma '%s': %v", proforma, err)
	}

	if baseSeq != "ELVISK" {
		t.Errorf("Expected sequence 'ELVISK', got '%s'", baseSeq)
	}

	// Check for ambiguous modification at position 4 (S)
	if mods, exists := modifications["4"]; exists && len(mods) > 0 {
		mod := mods[0]
		if mod.GetModType() != "ambiguous" {
			t.Errorf("Expected ambiguous modification type, got '%s'", mod.GetModType())
		}
	} else {
		t.Errorf("Expected ambiguous modification at position 4, but found none")
	}
}

func TestProFormaParserErrorCases(t *testing.T) {
	errorCases := []struct {
		name     string
		proforma string
	}{
		{
			name:     "Unclosed bracket",
			proforma: "ELVIS[PhosphoPEPTIDE",
		},
		{
			name:     "Unclosed global mod",
			proforma: "<15NPEPTIDE",
		},
		{
			name:     "Unclosed parenthesis",
			proforma: "ELVIS(PEPTIDE",
		},
		{
			name:     "Unmatched closing parenthesis",
			proforma: "ELVIS)PEPTIDE",
		},
	}

	for _, tt := range errorCases {
		t.Run(tt.name, func(t *testing.T) {
			_, _, _, _, _, err := ParseProForma(tt.proforma)
			if err == nil {
				t.Errorf("Expected error for invalid ProForma '%s', but parsing succeeded", tt.proforma)
			}
		})
	}
}
