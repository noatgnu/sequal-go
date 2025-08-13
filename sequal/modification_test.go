package sequal

import (
	"testing"
)

func TestModificationCreation(t *testing.T) {
	tests := []struct {
		name         string
		value        string
		modType      string
		expectedType string
		expectedValue string
	}{
		{
			name:         "Simple modification",
			value:        "Phospho",
			modType:      "static",
			expectedType: "static",
			expectedValue: "Phospho",
		},
		{
			name:         "Mass shift modification",
			value:        "+79.966",
			modType:      "static",
			expectedType: "static",
			expectedValue: "+79.966",
		},
		{
			name:         "Terminal modification",
			value:        "Acetyl",
			modType:      "terminal",
			expectedType: "terminal",
			expectedValue: "Acetyl",
		},
		{
			name:         "Labile modification",
			value:        "Glycan:Hex(1)HexNAc(2)",
			modType:      "labile",
			expectedType: "labile",
			expectedValue: "Hex(1)HexNAc(2)",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			mod := NewModification(tt.value, nil, nil, nil, tt.modType, false, 0, 0.0, false,
				nil, false, false, false, nil, false, false, nil, nil, nil, nil)

			if mod.GetModType() != tt.expectedType {
				t.Errorf("Expected mod type '%s', got '%s'", tt.expectedType, mod.GetModType())
			}

			if mod.GetValue() != tt.expectedValue {
				t.Errorf("Expected value '%s', got '%s'", tt.expectedValue, mod.GetValue())
			}
		})
	}
}

func TestModificationWithMass(t *testing.T) {
	mass := 79.966331
	mod := NewModification("Phospho", nil, nil, nil, "static", false, 0, mass, false,
		nil, false, false, false, nil, false, false, nil, nil, nil, nil)

	retrievedMass := mod.GetMass()
	if retrievedMass == nil {
		t.Errorf("Expected mass %f, got nil", mass)
	} else if *retrievedMass != mass {
		t.Errorf("Expected mass %f, got %f", mass, *retrievedMass)
	}
}

func TestModificationWithCrosslink(t *testing.T) {
	crosslinkID := "XL1"
	mod := NewModification("DSS", nil, nil, nil, "crosslink", false, 0, 138.068, false,
		&crosslinkID, false, false, false, nil, false, false, nil, nil, nil, nil)

	if mod.GetModType() != "crosslink" {
		t.Errorf("Expected mod type 'crosslink', got '%s'", mod.GetModType())
	}

	retrievedID := mod.GetCrosslinkID()
	if retrievedID == nil {
		t.Errorf("Expected crosslink ID '%s', got nil", crosslinkID)
	} else if *retrievedID != crosslinkID {
		t.Errorf("Expected crosslink ID '%s', got '%s'", crosslinkID, *retrievedID)
	}
}

func TestModificationWithAmbiguity(t *testing.T) {
	ambiguityGroup := "1"
	localizationScore := 0.75
	mod := NewModification("Phospho", nil, nil, nil, "ambiguous", false, 0, 0.0, false,
		nil, false, false, false, &ambiguityGroup, false, false, nil, nil, &localizationScore, nil)

	if mod.GetModType() != "ambiguous" {
		t.Errorf("Expected mod type 'ambiguous', got '%s'", mod.GetModType())
	}

	retrievedGroup := mod.GetAmbiguityGroup()
	if retrievedGroup == nil {
		t.Errorf("Expected ambiguity group '%s', got nil", ambiguityGroup)
	} else if *retrievedGroup != ambiguityGroup {
		t.Errorf("Expected ambiguity group '%s', got '%s'", ambiguityGroup, *retrievedGroup)
	}
}

func TestModificationToProforma(t *testing.T) {
	tests := []struct {
		name           string
		mod            *Modification
		expectedOutput string
	}{
		{
			name: "Simple modification",
			mod: NewModification("Phospho", nil, nil, nil, "static", false, 0, 0.0, false,
				nil, false, false, false, nil, false, false, nil, nil, nil, nil),
			expectedOutput: "Phospho",
		},
		{
			name: "Mass shift modification",
			mod: NewModification("+79.966", nil, nil, nil, "static", false, 0, 79.966, false,
				nil, false, false, false, nil, false, false, nil, nil, nil, nil),
			expectedOutput: "+79.966",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			output := tt.mod.ToProforma()
			if output != tt.expectedOutput {
				t.Errorf("Expected ProForma output '%s', got '%s'", tt.expectedOutput, output)
			}
		})
	}
}

func TestModificationEquality(t *testing.T) {
	mod1 := NewModification("Phospho", IntPtr(5), nil, nil, "static", false, 0, 79.966, false,
		nil, false, false, false, nil, false, false, nil, nil, nil, nil)

	mod2 := NewModification("Phospho", IntPtr(5), nil, nil, "static", false, 0, 79.966, false,
		nil, false, false, false, nil, false, false, nil, nil, nil, nil)

	mod3 := NewModification("Acetyl", IntPtr(5), nil, nil, "static", false, 0, 42.011, false,
		nil, false, false, false, nil, false, false, nil, nil, nil, nil)

	if !mod1.Equal(*mod2) {
		t.Errorf("Expected mod1 and mod2 to be equal")
	}

	if mod1.Equal(*mod3) {
		t.Errorf("Expected mod1 and mod3 to be different")
	}
}

func TestModificationToMap(t *testing.T) {
	mod := NewModification("Phospho", IntPtr(5), nil, StringPtr("Phosphorylation"), "static", false, 0, 79.966, false,
		nil, false, false, false, nil, false, false, nil, nil, nil, nil)

	modMap := mod.ToMap()

	if modMap["value"] != "Phospho" {
		t.Errorf("Expected value 'Phospho', got %v", modMap["value"])
	}

	if modMap["mod_type"] != "static" {
		t.Errorf("Expected mod_type 'static', got %v", modMap["mod_type"])
	}

	if modMap["full_name"] == nil {
		t.Errorf("Expected full_name to be present")
	} else if *modMap["full_name"].(*string) != "Phosphorylation" {
		t.Errorf("Expected full_name 'Phosphorylation', got %v", *modMap["full_name"].(*string))
	}
}

func TestModificationString(t *testing.T) {
	tests := []struct {
		name           string
		mod            *Modification
		expectedString string
	}{
		{
			name: "Simple modification",
			mod: NewModification("Phospho", nil, nil, nil, "static", false, 0, 0.0, false,
				nil, false, false, false, nil, false, false, nil, nil, nil, nil),
			expectedString: "Phospho",
		},
		{
			name: "Crosslink reference",
			mod: NewModification("", nil, nil, nil, "crosslink", false, 0, 0.0, false,
				StringPtr("XL1"), true, false, false, nil, false, false, nil, nil, nil, nil),
			expectedString: "#XL1",
		},
		{
			name: "Branch reference",
			mod: NewModification("", nil, nil, nil, "branch", false, 0, 0.0, false,
				nil, false, true, false, nil, false, false, nil, nil, nil, nil),
			expectedString: "#BRANCH",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result := tt.mod.String()
			if result != tt.expectedString {
				t.Errorf("Expected string '%s', got '%s'", tt.expectedString, result)
			}
		})
	}
}

func TestLabileModification(t *testing.T) {
	mod := NewModification("Glycan:Hex(1)HexNAc(2)", nil, nil, nil, "labile", true, 1, 0.0, false,
		nil, false, false, false, nil, false, false, nil, nil, nil, nil)

	if !mod.IsLabile() {
		t.Errorf("Expected modification to be labile")
	}

	if mod.GetLabileNumber() != 1 {
		t.Errorf("Expected labile number 1, got %d", mod.GetLabileNumber())
	}

	if mod.GetModType() != "labile" {
		t.Errorf("Expected mod type 'labile', got '%s'", mod.GetModType())
	}
}

func TestAllFilledModification(t *testing.T) {
	mod := NewModification("Carbamidomethyl", nil, nil, nil, "static", false, 0, 57.021, true,
		nil, false, false, false, nil, false, false, nil, nil, nil, nil)

	if !mod.IsAllFilled() {
		t.Errorf("Expected modification to be all-filled")
	}
}

func TestModificationWithRegex(t *testing.T) {
	regexPattern := "[ST]"
	mod := NewModification("Phospho", nil, &regexPattern, nil, "variable", false, 0, 79.966, false,
		nil, false, false, false, nil, false, false, nil, nil, nil, nil)

	regex := mod.GetRegex()
	if regex == nil {
		t.Errorf("Expected regex pattern to be set")
	} else {
		// Test the regex works
		testSequence := "PEPTIDES"
		positions := mod.FindPositions(testSequence)
		if len(positions) == 0 {
			t.Errorf("Expected to find positions for S and T residues")
		}

		// Should find S at position 6 and T at position 3 (0-indexed)
		foundS := false
		foundT := false
		for _, pos := range positions {
			if len(pos) >= 2 {
				if pos[0] == 7 { // S position
					foundS = true
				}
				if pos[0] == 3 { // T position  
					foundT = true
				}
			}
		}

		if !foundS {
			t.Errorf("Expected to find S residue at position 7")
		}
		if !foundT {
			t.Errorf("Expected to find T residue at position 3")
		}
	}
}

func TestModificationHash(t *testing.T) {
	mod1 := NewModification("Phospho", IntPtr(5), nil, nil, "static", false, 0, 79.966, false,
		nil, false, false, false, nil, false, false, nil, nil, nil, nil)

	mod2 := NewModification("Phospho", IntPtr(5), nil, nil, "static", false, 0, 79.966, false,
		nil, false, false, false, nil, false, false, nil, nil, nil, nil)

	hash1, err1 := mod1.Hash()
	hash2, err2 := mod2.Hash()

	if err1 != nil {
		t.Errorf("Failed to generate hash for mod1: %v", err1)
	}
	if err2 != nil {
		t.Errorf("Failed to generate hash for mod2: %v", err2)
	}

	if hash1 != hash2 {
		t.Errorf("Expected equal modifications to have same hash, got '%s' and '%s'", hash1, hash2)
	}

	// Different modification should have different hash
	mod3 := NewModification("Acetyl", IntPtr(5), nil, nil, "static", false, 0, 42.011, false,
		nil, false, false, false, nil, false, false, nil, nil, nil, nil)

	hash3, err3 := mod3.Hash()
	if err3 != nil {
		t.Errorf("Failed to generate hash for mod3: %v", err3)
	}

	if hash1 == hash3 {
		t.Errorf("Expected different modifications to have different hashes")
	}
}