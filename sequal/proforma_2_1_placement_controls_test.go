package sequal

import (
	"testing"
)

func TestPlacementControls_PositionConstraint(t *testing.T) {
	tests := []struct {
		name           string
		proformaString string
		expectedPos    []string
		shouldParse    bool
	}{
		{
			name:           "Single position constraint",
			proformaString: "<[Oxidation|Position:M]@M>PEPTIDE",
			expectedPos:    []string{"M"},
			shouldParse:    true,
		},
		{
			name:           "Multiple position constraints",
			proformaString: "<[Phospho|Position:S,T,Y]@S,T,Y>PEPTIDES",
			expectedPos:    []string{"S", "T", "Y"},
			shouldParse:    true,
		},
		{
			name:           "Position constraint with specific amino acid",
			proformaString: "<[Carbamidomethyl|Position:C]@C>PEPTCDE",
			expectedPos:    []string{"C"},
			shouldParse:    true,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			seq, err := FromProforma(tt.proformaString)
			if tt.shouldParse {
				if err != nil {
					t.Errorf("Failed to parse %s: %v", tt.proformaString, err)
					return
				}

				if len(seq.GetGlobalMods()) == 0 {
					t.Errorf("Expected global modification for %s", tt.proformaString)
					return
				}

				mod := seq.GetGlobalMods()[0]
				posConstraint := mod.GetPositionConstraint()

				if posConstraint == nil {
					t.Errorf("Expected position constraint for %s", tt.proformaString)
					return
				}

				if len(posConstraint) != len(tt.expectedPos) {
					t.Errorf("Expected %d positions, got %d for %s", len(tt.expectedPos), len(posConstraint), tt.proformaString)
					return
				}

				for i, expectedP := range tt.expectedPos {
					if posConstraint[i] != expectedP {
						t.Errorf("Expected position[%d]='%s', got '%s' for %s", i, expectedP, posConstraint[i], tt.proformaString)
					}
				}
			} else {
				if err == nil {
					t.Errorf("Expected error parsing %s", tt.proformaString)
				}
			}
		})
	}
}

func TestPlacementControls_LimitPerPosition(t *testing.T) {
	tests := []struct {
		name           string
		proformaString string
		expectedLimit  int
		shouldParse    bool
	}{
		{
			name:           "Limit of 1",
			proformaString: "<[Oxidation|Limit:1]@M>MMMM",
			expectedLimit:  1,
			shouldParse:    true,
		},
		{
			name:           "Limit of 2",
			proformaString: "<[Phospho|Limit:2]@S,T,Y>STYSTY",
			expectedLimit:  2,
			shouldParse:    true,
		},
		{
			name:           "Limit with position constraint",
			proformaString: "<[Oxidation|Position:M|Limit:1]@M>MMMM",
			expectedLimit:  1,
			shouldParse:    true,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			seq, err := FromProforma(tt.proformaString)
			if tt.shouldParse {
				if err != nil {
					t.Errorf("Failed to parse %s: %v", tt.proformaString, err)
					return
				}

				if len(seq.GetGlobalMods()) == 0 {
					t.Errorf("Expected global modification for %s", tt.proformaString)
					return
				}

				mod := seq.GetGlobalMods()[0]
				limit := mod.GetLimitPerPosition()

				if limit == nil {
					t.Errorf("Expected limit per position for %s", tt.proformaString)
					return
				}

				if *limit != tt.expectedLimit {
					t.Errorf("Expected limit=%d, got %d for %s", tt.expectedLimit, *limit, tt.proformaString)
				}
			} else {
				if err == nil {
					t.Errorf("Expected error parsing %s", tt.proformaString)
				}
			}
		})
	}
}

func TestPlacementControls_CoMKP(t *testing.T) {
	tests := []struct {
		name            string
		proformaString  string
		shouldHaveCoMKP bool
		shouldParse     bool
	}{
		{
			name:            "CoMKP tag present",
			proformaString:  "<[Oxidation|CoMKP]@M>PEPTIDE",
			shouldHaveCoMKP: true,
			shouldParse:     true,
		},
		{
			name:            "Full form ColocaliseModificationsOfKnownPosition",
			proformaString:  "<[Oxidation|ColocaliseModificationsOfKnownPosition]@M>PEPTIDE",
			shouldHaveCoMKP: true,
			shouldParse:     true,
		},
		{
			name:            "CoMKP with other tags",
			proformaString:  "<[Phospho|Position:S,T|CoMKP]@S,T>STYSTY",
			shouldHaveCoMKP: true,
			shouldParse:     true,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			seq, err := FromProforma(tt.proformaString)
			if tt.shouldParse {
				if err != nil {
					t.Errorf("Failed to parse %s: %v", tt.proformaString, err)
					return
				}

				if len(seq.GetGlobalMods()) == 0 {
					t.Errorf("Expected global modification for %s", tt.proformaString)
					return
				}

				mod := seq.GetGlobalMods()[0]
				hasCoMKP := mod.GetColocalizeKnown()

				if hasCoMKP != tt.shouldHaveCoMKP {
					t.Errorf("Expected CoMKP=%v, got %v for %s", tt.shouldHaveCoMKP, hasCoMKP, tt.proformaString)
				}
			} else {
				if err == nil {
					t.Errorf("Expected error parsing %s", tt.proformaString)
				}
			}
		})
	}
}

func TestPlacementControls_CoMUP(t *testing.T) {
	tests := []struct {
		name            string
		proformaString  string
		shouldHaveCoMUP bool
		shouldParse     bool
	}{
		{
			name:            "CoMUP tag present",
			proformaString:  "<[Dioxidation|CoMUP]@M>PEPTIDE",
			shouldHaveCoMUP: true,
			shouldParse:     true,
		},
		{
			name:            "Full form ColocaliseModificationsOfUnknownPosition",
			proformaString:  "<[Dioxidation|ColocaliseModificationsOfUnknownPosition]@M>PEPTIDE",
			shouldHaveCoMUP: true,
			shouldParse:     true,
		},
		{
			name:            "CoMUP with other tags",
			proformaString:  "<[Oxidation|Position:M|CoMUP]@M>PEPTIDE",
			shouldHaveCoMUP: true,
			shouldParse:     true,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			seq, err := FromProforma(tt.proformaString)
			if tt.shouldParse {
				if err != nil {
					t.Errorf("Failed to parse %s: %v", tt.proformaString, err)
					return
				}

				if len(seq.GetGlobalMods()) == 0 {
					t.Errorf("Expected global modification for %s", tt.proformaString)
					return
				}

				mod := seq.GetGlobalMods()[0]
				hasCoMUP := mod.GetColocalizeUnknown()

				if hasCoMUP != tt.shouldHaveCoMUP {
					t.Errorf("Expected CoMUP=%v, got %v for %s", tt.shouldHaveCoMUP, hasCoMUP, tt.proformaString)
				}
			} else {
				if err == nil {
					t.Errorf("Expected error parsing %s", tt.proformaString)
				}
			}
		})
	}
}

func TestPlacementControls_CombinedTags(t *testing.T) {
	tests := []struct {
		name           string
		proformaString string
		description    string
		shouldParse    bool
	}{
		{
			name:           "All placement controls together",
			proformaString: "<[Phospho|Position:S,T,Y|Limit:2|CoMKP]@S,T,Y>STYSTY",
			description:    "Position, Limit, and CoMKP together",
			shouldParse:    true,
		},
		{
			name:           "Position and Limit",
			proformaString: "<[Oxidation|Position:M|Limit:1]@M>MMMM",
			description:    "Position and Limit together",
			shouldParse:    true,
		},
		{
			name:           "Position and CoMUP",
			proformaString: "<[Acetyl|Position:K|CoMUP]@K>KKKK",
			description:    "Position and CoMUP together",
			shouldParse:    true,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			seq, err := FromProforma(tt.proformaString)
			if tt.shouldParse {
				if err != nil {
					t.Errorf("Failed to parse %s: %v", tt.proformaString, err)
					return
				}

				if len(seq.GetGlobalMods()) == 0 {
					t.Errorf("Expected global modification for %s", tt.proformaString)
					return
				}

				mod := seq.GetGlobalMods()[0]

				if tt.name == "All placement controls together" {
					if mod.GetPositionConstraint() == nil || len(mod.GetPositionConstraint()) != 3 {
						t.Errorf("Expected 3 position constraints")
					}
					if mod.GetLimitPerPosition() == nil || *mod.GetLimitPerPosition() != 2 {
						t.Errorf("Expected limit=2")
					}
					if !mod.GetColocalizeKnown() {
						t.Errorf("Expected CoMKP=true")
					}
				}
			} else {
				if err == nil {
					t.Errorf("Expected error parsing %s", tt.proformaString)
				}
			}
		})
	}
}

func TestPlacementControls_Serialization(t *testing.T) {
	tests := []struct {
		name           string
		posConstraint  []string
		limitPerPos    *int
		coMKP          bool
		coMUP          bool
		expectedOutput string
	}{
		{
			name:           "Position constraint only",
			posConstraint:  []string{"M"},
			limitPerPos:    nil,
			coMKP:          false,
			coMUP:          false,
			expectedOutput: "Phospho|Position:M",
		},
		{
			name:           "Multiple positions",
			posConstraint:  []string{"S", "T", "Y"},
			limitPerPos:    nil,
			coMKP:          false,
			coMUP:          false,
			expectedOutput: "Phospho|Position:S,T,Y",
		},
		{
			name:           "Limit only",
			posConstraint:  nil,
			limitPerPos:    IntPtr(2),
			coMKP:          false,
			coMUP:          false,
			expectedOutput: "Phospho|Limit:2",
		},
		{
			name:           "CoMKP only",
			posConstraint:  nil,
			limitPerPos:    nil,
			coMKP:          true,
			coMUP:          false,
			expectedOutput: "Phospho|CoMKP",
		},
		{
			name:           "CoMUP only",
			posConstraint:  nil,
			limitPerPos:    nil,
			coMKP:          false,
			coMUP:          true,
			expectedOutput: "Phospho|CoMUP",
		},
		{
			name:           "All combined",
			posConstraint:  []string{"S", "T"},
			limitPerPos:    IntPtr(1),
			coMKP:          true,
			coMUP:          false,
			expectedOutput: "Phospho|Position:S,T|Limit:1|CoMKP",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			mod := NewModification("Phospho", nil, nil, nil, "static", false, 0, 79.966, false,
				nil, false, false, false, nil, false, false, nil, nil, nil, nil,
				tt.posConstraint, tt.limitPerPos, tt.coMKP, tt.coMUP, false)

			output := mod.ToProforma()
			if output != tt.expectedOutput {
				t.Errorf("Expected '%s', got '%s'", tt.expectedOutput, output)
			}
		})
	}
}

func TestPlacementControls_RoundTrip(t *testing.T) {
	tests := []struct {
		name           string
		proformaString string
	}{
		{
			name:           "Position constraint roundtrip",
			proformaString: "<[Oxidation|Position:M]@M>MTPEPTIDE",
		},
		{
			name:           "Limit roundtrip",
			proformaString: "<[Phospho|Limit:2]@S,T,Y>PEPTIDES",
		},
		{
			name:           "CoMKP roundtrip",
			proformaString: "<[Oxidation|CoMKP]@M>PEPTIDE",
		},
		{
			name:           "Combined roundtrip",
			proformaString: "<[Phospho|Position:S,T|Limit:1|CoMKP]@S,T>PEPTIDES",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			seq, err := FromProforma(tt.proformaString)
			if err != nil {
				t.Errorf("Failed to parse %s: %v", tt.proformaString, err)
				return
			}

			output := seq.ToProforma()

			seq2, err := FromProforma(output)
			if err != nil {
				t.Errorf("Failed to re-parse output %s: %v", output, err)
				return
			}

			if len(seq.GetGlobalMods()) != len(seq2.GetGlobalMods()) {
				t.Errorf("Global modification count mismatch after roundtrip")
				return
			}

			if len(seq.GetGlobalMods()) > 0 {
				mod1 := seq.GetGlobalMods()[0]
				mod2 := seq2.GetGlobalMods()[0]

				pos1 := mod1.GetPositionConstraint()
				pos2 := mod2.GetPositionConstraint()
				if (pos1 == nil) != (pos2 == nil) {
					t.Errorf("Position constraint presence mismatch after roundtrip")
				}

				limit1 := mod1.GetLimitPerPosition()
				limit2 := mod2.GetLimitPerPosition()
				if (limit1 == nil) != (limit2 == nil) {
					t.Errorf("Limit presence mismatch after roundtrip")
				}

				if mod1.GetColocalizeKnown() != mod2.GetColocalizeKnown() {
					t.Errorf("CoMKP mismatch after roundtrip")
				}

				if mod1.GetColocalizeUnknown() != mod2.GetColocalizeUnknown() {
					t.Errorf("CoMUP mismatch after roundtrip")
				}
			}
		})
	}
}
