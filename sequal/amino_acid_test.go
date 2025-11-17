package sequal

import (
	"testing"
)

func TestAminoAcidCreation(t *testing.T) {
	tests := []struct {
		name          string
		value         string
		position      *int
		mass          *float64
		expectedValue string
		expectedMass  float64
		shouldError   bool
	}{
		{
			name:          "Alanine",
			value:         "A",
			position:      IntPtr(0),
			mass:          nil,
			expectedValue: "A",
			expectedMass:  71.037114,
			shouldError:   false,
		},
		{
			name:          "Serine",
			value:         "S",
			position:      IntPtr(5),
			mass:          nil,
			expectedValue: "S",
			expectedMass:  87.032028,
			shouldError:   false,
		},
		{
			name:          "Custom mass",
			value:         "X",
			position:      IntPtr(1),
			mass:          Float64Ptr(123.456),
			expectedValue: "X",
			expectedMass:  123.456,
			shouldError:   false,
		},
		{
			name:          "Unknown without mass",
			value:         "Z",
			position:      IntPtr(1),
			mass:          nil,
			expectedValue: "Z",
			expectedMass:  0,
			shouldError:   true,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			aa, err := NewAminoAcid(tt.value, tt.position, tt.mass)

			if tt.shouldError {
				if err == nil {
					t.Errorf("Expected error for unknown amino acid without mass, but got none")
				}
				return
			}

			if err != nil {
				t.Fatalf("Unexpected error creating amino acid: %v", err)
			}

			if aa.GetValue() != tt.expectedValue {
				t.Errorf("Expected value '%s', got '%s'", tt.expectedValue, aa.GetValue())
			}

			retrievedMass := aa.GetMass()
			if retrievedMass == nil {
				t.Errorf("Expected mass %f, got nil", tt.expectedMass)
			} else if *retrievedMass != tt.expectedMass {
				t.Errorf("Expected mass %f, got %f", tt.expectedMass, *retrievedMass)
			}

			if aa.GetPosition() == nil {
				if tt.position != nil {
					t.Errorf("Expected position %d, got nil", *tt.position)
				}
			} else if tt.position == nil {
				t.Errorf("Expected nil position, got %d", *aa.GetPosition())
			} else if *aa.GetPosition() != *tt.position {
				t.Errorf("Expected position %d, got %d", *tt.position, *aa.GetPosition())
			}
		})
	}
}

func TestAminoAcidModifications(t *testing.T) {
	aa, err := NewAminoAcid("S", IntPtr(5), nil)
	if err != nil {
		t.Fatalf("Failed to create amino acid: %v", err)
	}

	// Initially no modifications
	mods := aa.GetMods()
	if len(mods) != 0 {
		t.Errorf("Expected no modifications initially, got %d", len(mods))
	}

	// Add a modification
	phosphoMod := NewModification("Phospho", nil, nil, nil, "static", false, 0, 79.966, false,
		nil, false, false, false, nil, false, false, nil, nil, nil, nil,
		nil, nil, false, false, false)

	aa.AddModification(phosphoMod)

	// Check modification was added
	mods = aa.GetMods()
	if len(mods) != 1 {
		t.Errorf("Expected 1 modification, got %d", len(mods))
	}

	if mods[0].GetValue() != "Phospho" {
		t.Errorf("Expected modification 'Phospho', got '%s'", mods[0].GetValue())
	}

	// Check HasModification
	if !aa.HasModification("Phospho") {
		t.Errorf("Expected amino acid to have 'Phospho' modification")
	}

	if aa.HasModification("Acetyl") {
		t.Errorf("Expected amino acid NOT to have 'Acetyl' modification")
	}

	// Test HasModification with Modification object
	if !aa.HasModification(phosphoMod) {
		t.Errorf("Expected amino acid to have the phospho modification object")
	}

	// Add another modification
	acetylMod := NewModification("Acetyl", nil, nil, nil, "static", false, 0, 42.011, false,
		nil, false, false, false, nil, false, false, nil, nil, nil, nil,
		nil, nil, false, false, false)

	aa.AddModification(acetylMod)

	mods = aa.GetMods()
	if len(mods) != 2 {
		t.Errorf("Expected 2 modifications, got %d", len(mods))
	}
}

func TestAminoAcidRemoveModifications(t *testing.T) {
	aa, err := NewAminoAcid("S", IntPtr(5), nil)
	if err != nil {
		t.Fatalf("Failed to create amino acid: %v", err)
	}

	// Add modifications
	phosphoMod := NewModification("Phospho", nil, nil, nil, "static", false, 0, 79.966, false,
		nil, false, false, false, nil, false, false, nil, nil, nil, nil,
		nil, nil, false, false, false)
	acetylMod := NewModification("Acetyl", nil, nil, nil, "static", false, 0, 42.011, false,
		nil, false, false, false, nil, false, false, nil, nil, nil, nil,
		nil, nil, false, false, false)

	aa.AddModification(phosphoMod)
	aa.AddModification(acetylMod)

	if len(aa.GetMods()) != 2 {
		t.Errorf("Expected 2 modifications initially, got %d", len(aa.GetMods()))
	}

	// Remove by string
	removed := aa.RemoveModification("Phospho")
	if !removed {
		t.Errorf("Expected to successfully remove 'Phospho' modification")
	}

	if len(aa.GetMods()) != 1 {
		t.Errorf("Expected 1 modification after removal, got %d", len(aa.GetMods()))
	}

	if aa.HasModification("Phospho") {
		t.Errorf("Expected 'Phospho' modification to be removed")
	}

	// Remove by object
	removed = aa.RemoveModification(acetylMod)
	if !removed {
		t.Errorf("Expected to successfully remove acetyl modification object")
	}

	if len(aa.GetMods()) != 0 {
		t.Errorf("Expected 0 modifications after removal, got %d", len(aa.GetMods()))
	}

	// Try to remove non-existent modification
	removed = aa.RemoveModification("NonExistent")
	if removed {
		t.Errorf("Expected NOT to remove non-existent modification")
	}
}

func TestAminoAcidTotalMass(t *testing.T) {
	aa, err := NewAminoAcid("S", IntPtr(5), nil)
	if err != nil {
		t.Fatalf("Failed to create amino acid: %v", err)
	}

	// Initial mass should be just the amino acid mass
	expectedInitialMass := AAMass["S"]
	totalMass := aa.GetTotalMass()
	if totalMass != expectedInitialMass {
		t.Errorf("Expected initial mass %f, got %f", expectedInitialMass, totalMass)
	}

	// Add modifications with mass
	phosphoMod := NewModification("Phospho", nil, nil, nil, "static", false, 0, 79.966, false,
		nil, false, false, false, nil, false, false, nil, nil, nil, nil,
		nil, nil, false, false, false)
	aa.AddModification(phosphoMod)

	expectedTotalMass := expectedInitialMass + 79.966
	totalMass = aa.GetTotalMass()
	if totalMass != expectedTotalMass {
		t.Errorf("Expected total mass %f, got %f", expectedTotalMass, totalMass)
	}

	// Add another modification
	acetylMod := NewModification("Acetyl", nil, nil, nil, "static", false, 0, 42.011, false,
		nil, false, false, false, nil, false, false, nil, nil, nil, nil,
		nil, nil, false, false, false)
	aa.AddModification(acetylMod)

	expectedTotalMass += 42.011
	totalMass = aa.GetTotalMass()
	if totalMass != expectedTotalMass {
		t.Errorf("Expected total mass %f, got %f", expectedTotalMass, totalMass)
	}
}

func TestAminoAcidToMap(t *testing.T) {
	aa, err := NewAminoAcid("S", IntPtr(5), nil)
	if err != nil {
		t.Fatalf("Failed to create amino acid: %v", err)
	}

	// Add a modification
	phosphoMod := NewModification("Phospho", nil, nil, nil, "static", false, 0, 79.966, false,
		nil, false, false, false, nil, false, false, nil, nil, nil, nil,
		nil, nil, false, false, false)
	aa.AddModification(phosphoMod)

	aaMap := aa.ToMap()

	if aaMap["value"] != "S" {
		t.Errorf("Expected value 'S', got %v", aaMap["value"])
	}

	if aaMap["position"] == nil {
		t.Errorf("Expected position to be present")
	} else if *aaMap["position"].(*int) != 5 {
		t.Errorf("Expected position 5, got %v", *aaMap["position"].(*int))
	}

	expectedTotalMass := AAMass["S"] + 79.966
	if aaMap["total_mass"] != expectedTotalMass {
		t.Errorf("Expected total mass %f, got %v", expectedTotalMass, aaMap["total_mass"])
	}

	// Check modifications are included
	if aaMap["mods"] == nil {
		t.Errorf("Expected mods to be present in map")
	} else {
		modsSlice := aaMap["mods"].([]map[string]interface{})
		if len(modsSlice) != 1 {
			t.Errorf("Expected 1 modification in map, got %d", len(modsSlice))
		}
	}
}

func TestAminoAcidEquality(t *testing.T) {
	aa1, err1 := NewAminoAcid("S", IntPtr(5), nil)
	aa2, err2 := NewAminoAcid("S", IntPtr(5), nil)
	aa3, err3 := NewAminoAcid("T", IntPtr(5), nil)

	if err1 != nil || err2 != nil || err3 != nil {
		t.Fatalf("Failed to create amino acids: %v, %v, %v", err1, err2, err3)
	}

	// Same amino acids should be equal
	if !aa1.Equal(aa2) {
		t.Errorf("Expected aa1 and aa2 to be equal")
	}

	// Different amino acids should not be equal
	if aa1.Equal(aa3) {
		t.Errorf("Expected aa1 and aa3 to be different")
	}

	// Add modifications and test equality
	phosphoMod1 := NewModification("Phospho", nil, nil, nil, "static", false, 0, 79.966, false,
		nil, false, false, false, nil, false, false, nil, nil, nil, nil,
		nil, nil, false, false, false)
	phosphoMod2 := NewModification("Phospho", nil, nil, nil, "static", false, 0, 79.966, false,
		nil, false, false, false, nil, false, false, nil, nil, nil, nil,
		nil, nil, false, false, false)

	aa1.AddModification(phosphoMod1)
	aa2.AddModification(phosphoMod2)

	if !aa1.Equal(aa2) {
		t.Errorf("Expected aa1 and aa2 to be equal even with modifications")
	}

	// Add different modification to one
	acetylMod := NewModification("Acetyl", nil, nil, nil, "static", false, 0, 42.011, false,
		nil, false, false, false, nil, false, false, nil, nil, nil, nil,
		nil, nil, false, false, false)
	aa2.AddModification(acetylMod)

	if aa1.Equal(aa2) {
		t.Errorf("Expected aa1 and aa2 to be different with different modifications")
	}
}

func TestAminoAcidString(t *testing.T) {
	aa, err := NewAminoAcid("S", IntPtr(5), nil)
	if err != nil {
		t.Fatalf("Failed to create amino acid: %v", err)
	}

	// Without modifications
	str := aa.String()
	if str != "S" {
		t.Errorf("Expected string 'S', got '%s'", str)
	}

	// With modifications
	phosphoMod := NewModification("Phospho", nil, nil, nil, "static", false, 0, 79.966, false,
		nil, false, false, false, nil, false, false, nil, nil, nil, nil,
		nil, nil, false, false, false)
	aa.AddModification(phosphoMod)

	str = aa.String()
	expected := "S[Phospho]"
	if str != expected {
		t.Errorf("Expected string '%s', got '%s'", expected, str)
	}

	// With multiple modifications
	acetylMod := NewModification("Acetyl", nil, nil, nil, "static", false, 0, 42.011, false,
		nil, false, false, false, nil, false, false, nil, nil, nil, nil,
		nil, nil, false, false, false)
	aa.AddModification(acetylMod)

	str = aa.String()
	// Should contain both modifications
	if !containsString(str, "[Phospho]") {
		t.Errorf("Expected string to contain '[Phospho]', got '%s'", str)
	}
	if !containsString(str, "[Acetyl]") {
		t.Errorf("Expected string to contain '[Acetyl]', got '%s'", str)
	}
}

func TestAminoAcidHashCode(t *testing.T) {
	aa1, err1 := NewAminoAcid("S", IntPtr(5), nil)
	aa2, err2 := NewAminoAcid("S", IntPtr(5), nil)
	aa3, err3 := NewAminoAcid("T", IntPtr(5), nil)

	if err1 != nil || err2 != nil || err3 != nil {
		t.Fatalf("Failed to create amino acids: %v, %v, %v", err1, err2, err3)
	}

	hash1, err1 := aa1.Hash()
	hash2, err2 := aa2.Hash()
	hash3, err3 := aa3.Hash()

	if err1 != nil || err2 != nil || err3 != nil {
		t.Fatalf("Failed to generate hashes: %v, %v, %v", err1, err2, err3)
	}

	// Same amino acids should have same hash
	if hash1 != hash2 {
		t.Errorf("Expected aa1 and aa2 to have same hash, got '%s' and '%s'", hash1, hash2)
	}

	// Different amino acids should have different hash
	if hash1 == hash3 {
		t.Errorf("Expected aa1 and aa3 to have different hashes")
	}
}

func TestAminoAcidDebugString(t *testing.T) {
	aa, err := NewAminoAcid("S", IntPtr(5), nil)
	if err != nil {
		t.Fatalf("Failed to create amino acid: %v", err)
	}

	debugStr := aa.ToDebugString()

	// Should contain value, position, and modification info
	if !containsString(debugStr, "S") {
		t.Errorf("Expected debug string to contain 'S', got '%s'", debugStr)
	}
	if !containsString(debugStr, "position=5") {
		t.Errorf("Expected debug string to contain 'position=5', got '%s'", debugStr)
	}
	if !containsString(debugStr, "AminoAcid") {
		t.Errorf("Expected debug string to contain 'AminoAcid', got '%s'", debugStr)
	}
}

func TestSetModification(t *testing.T) {
	// Test the legacy setModification method
	aa, err := NewAminoAcid("S", IntPtr(5), nil)
	if err != nil {
		t.Fatalf("Failed to create amino acid: %v", err)
	}

	phosphoMod := NewModification("Phospho", nil, nil, nil, "static", false, 0, 79.966, false,
		nil, false, false, false, nil, false, false, nil, nil, nil, nil,
		nil, nil, false, false, false)

	aa.SetModification(phosphoMod)

	mods := aa.GetMods()
	if len(mods) != 1 {
		t.Errorf("Expected 1 modification, got %d", len(mods))
	}

	if mods[0].GetValue() != "Phospho" {
		t.Errorf("Expected modification 'Phospho', got '%s'", mods[0].GetValue())
	}
}

// Helper function to check if a string contains a substring
func containsString(s, substr string) bool {
	for i := 0; i <= len(s)-len(substr); i++ {
		if s[i:i+len(substr)] == substr {
			return true
		}
	}
	return false
}
