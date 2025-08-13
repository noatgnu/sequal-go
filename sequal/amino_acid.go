package sequal

import (
	"crypto/sha256"
	"encoding/json"
	"fmt"
	"strings"
)

// AminoAcid represents an amino acid residue with its position, modifications, and mass properties.
// It embeds BaseBlock to inherit common functionality and maintains a list of modifications
// that can be applied to the residue.
type AminoAcid struct {
	BaseBlock
	mods []*Modification
}

// NewAminoAcid creates a new AminoAcid instance with the specified value, position, and optional mass.
// If the amino acid is not recognized and no mass is provided, it returns an error.
// The mass parameter overrides the default mass for known amino acids.
func NewAminoAcid(value string, position *int, mass *float64) (*AminoAcid, error) {
	if _, exists := AAMass[value]; !exists && mass == nil {
		return nil, fmt.Errorf("unknown amino acid '%s' and no mass provided", value)
	}

	var inferredMass float64
	if mass != nil {
		inferredMass = *mass
	} else {
		inferredMass = AAMass[value]
	}

	baseBlock := NewBaseBlock(value, position, false, &inferredMass)

	return &AminoAcid{
		BaseBlock: baseBlock,
		mods:      []*Modification{},
	}, nil
}

// GetMods returns a copy of the modifications list to prevent external mutation.
func (aa *AminoAcid) GetMods() []*Modification {
	modsCopy := make([]*Modification, len(aa.mods))
	copy(modsCopy, aa.mods)
	return modsCopy
}

// AddModification appends a modification to this amino acid's modification list.
func (aa *AminoAcid) AddModification(mod *Modification) {
	aa.mods = append(aa.mods, mod)
}

// SetModification adds a modification to this amino acid.
// This is a legacy method that calls AddModification for backward compatibility.
func (aa *AminoAcid) SetModification(mod *Modification) {
	aa.AddModification(mod)
}

// RemoveModification removes a modification from this amino acid.
// The mod parameter can be either a string (modification value) or a *Modification instance.
// Returns true if a modification was removed, false otherwise.
func (aa *AminoAcid) RemoveModification(mod interface{}) bool {
	switch m := mod.(type) {
	case string:
		for i, existingMod := range aa.mods {
			if existingMod.GetValue() == m {
				aa.mods = append(aa.mods[:i], aa.mods[i+1:]...)
				return true
			}
		}
	case *Modification:
		for i, existingMod := range aa.mods {
			if existingMod == m {
				aa.mods = append(aa.mods[:i], aa.mods[i+1:]...)
				return true
			}
		}
	}
	return false
}

// HasModification checks if this amino acid has a specific modification.
// The mod parameter can be either a string (modification value) or a *Modification instance.
func (aa *AminoAcid) HasModification(mod interface{}) bool {
	switch m := mod.(type) {
	case string:
		for _, existingMod := range aa.mods {
			if existingMod.GetValue() == m {
				return true
			}
		}
	case *Modification:
		for _, existingMod := range aa.mods {
			if existingMod == m {
				return true
			}
		}
	}
	return false
}

// GetTotalMass calculates the total mass of the amino acid including all modifications.
// Returns 0 if the base mass is not set.
func (aa *AminoAcid) GetTotalMass() float64 {
	mass := aa.GetMass()
	if mass == nil {
		return 0
	}

	total := *mass
	for _, mod := range aa.mods {
		if modMass := mod.GetMass(); modMass != nil {
			total += *modMass
		}
	}
	return total
}

// ToMap converts the amino acid to a map representation suitable for serialization.
// Includes base block properties, modifications, and total mass.
func (aa *AminoAcid) ToMap() map[string]interface{} {
	result := aa.BaseBlock.ToMap()

	mods := make([]map[string]interface{}, len(aa.mods))
	for i, mod := range aa.mods {
		mods[i] = mod.ToMap()
	}

	result["mods"] = mods
	result["total_mass"] = aa.GetTotalMass()

	return result
}

// Equal compares two amino acids for equality including their modifications.
// Returns true if both amino acids have the same value, position, mass, and modifications.
func (aa *AminoAcid) Equal(other BaseBlock) bool {
	otherAA, ok := other.(*AminoAcid)
	if !ok {
		return false
	}
	
	if aa.GetValue() != otherAA.GetValue() {
		return false
	}
	
	if aa.GetPosition() == nil && otherAA.GetPosition() == nil {
	} else if aa.GetPosition() != nil && otherAA.GetPosition() != nil {
		if *aa.GetPosition() != *otherAA.GetPosition() {
			return false
		}
	} else {
		return false
	}
	
	if aa.GetMass() == nil && otherAA.GetMass() == nil {
	} else if aa.GetMass() != nil && otherAA.GetMass() != nil {
		if *aa.GetMass() != *otherAA.GetMass() {
			return false
		}
	} else {
		return false
	}
	
	if len(aa.mods) != len(otherAA.mods) {
		return false
	}
	
	for i, mod := range aa.mods {
		if !mod.Equal(*otherAA.mods[i]) {
			return false
		}
	}
	
	return true
}

// Hash generates a SHA-256 hash for the amino acid including all modifications.
// The hash is computed from the JSON representation of the amino acid's map.
func (aa *AminoAcid) Hash() (string, error) {
	aaMap := aa.ToMap()
	jsonData, err := json.Marshal(aaMap)
	if err != nil {
		return "", err
	}

	modHash := sha256.Sum256(jsonData)

	return fmt.Sprintf("%x", modHash), nil
}

// ToDebugString returns a detailed string representation suitable for debugging.
// Includes the amino acid value, position, and all modifications.
func (aa *AminoAcid) ToDebugString() string {
	modStr := ""
	for i, mod := range aa.mods {
		if i > 0 {
			modStr += ", "
		}
		modStr += mod.String()
	}
	return fmt.Sprintf("AminoAcid(value='%s', position=%d, mods=[%s])", 
		aa.GetValue(), IntValue(aa.GetPosition()), modStr)
}

// String returns a string representation of the amino acid with its modifications
// in ProForma-like notation (e.g., "S[Phospho]").
func (aa *AminoAcid) String() string {
	var sb strings.Builder
	sb.WriteString(aa.GetValue())

	for _, mod := range aa.mods {
		sb.WriteString("[")
		sb.WriteString(mod.GetValue())
		sb.WriteString("]")
	}

	return sb.String()
}
