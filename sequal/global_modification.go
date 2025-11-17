package sequal

import (
	"fmt"
	"strings"
)

// GlobalModification represents a global modification that applies to specified residues
type GlobalModification struct {
	Modification
	targetResidues []string
	globalModType  string
}

// NewGlobalModification creates a new GlobalModification instance
func NewGlobalModification(value string, targetResidues []string, modType string, positionConstraint []string, limitPerPosition *int, colocalizeKnown bool, colocalizeUnknown bool) *GlobalModification {
	if modType != "isotope" && modType != "fixed" {
		panic("Global modification type must be 'isotope' or 'fixed'")
	}

	mod := NewModification(
		value,
		nil, // position
		nil, // regexPattern
		nil, // fullName
		"global",
		false,              // labile
		0,                  // labilNumber
		0,                  // mass
		false,              // allFilled
		nil,                // crosslinkID
		false,              // isCrosslinkRef
		false,              // isBranchRef
		false,              // isBranch
		nil,                // ambiguityGroup
		false,              // isAmbiguityRef
		false,              // inRange
		nil,                // rangeStart
		nil,                // rangeEnd
		nil,                // localizationScore
		nil,                // modValue
		positionConstraint, // positionConstraint (ProForma 2.1)
		limitPerPosition,   // limitPerPosition (ProForma 2.1)
		colocalizeKnown,    // colocalizeKnown (ProForma 2.1)
		colocalizeUnknown,  // colocalizeUnknown (ProForma 2.1)
		false,              // isIonType (ProForma 2.1)
	)

	return &GlobalModification{
		Modification:   *mod,
		targetResidues: targetResidues,
		globalModType:  modType,
	}
}

// GetTargetResidues returns the target residue types
func (gm *GlobalModification) GetTargetResidues() []string {
	return gm.targetResidues
}

// GetGlobalModType returns the global modification type
func (gm *GlobalModification) GetGlobalModType() string {
	return gm.globalModType
}

// ToProforma converts the modification to ProForma notation
func (gm *GlobalModification) ToProforma() string {
	if gm.globalModType == "isotope" {
		return fmt.Sprintf("<%s>", gm.Modification.ToProforma())
	} else {
		modValue := gm.Modification.ToProforma()
		var modStr string
		// Check if the modification is a mass (starts with + or -)
		// If not, it's a named modification and needs brackets
		if strings.HasPrefix(modValue, "+") || strings.HasPrefix(modValue, "-") {
			modStr = modValue
		} else {
			// Named modification - wrap in brackets
			modStr = fmt.Sprintf("[%s]", modValue)
		}
		targets := strings.Join(gm.targetResidues, ",")
		return fmt.Sprintf("<%s@%s>", modStr, targets)
	}
}

// String returns a string representation of the global modification
func (gm *GlobalModification) String() string {
	return gm.ToProforma()
}

// Repr returns a detailed string representation for debugging
func (gm *GlobalModification) Repr() string {
	base := fmt.Sprintf("GlobalModification(value='%s'", gm.GetValue())

	if gm.GetSource() != nil {
		base += fmt.Sprintf(", source='%s'", *gm.GetSource())
	}

	if gm.targetResidues != nil {
		base += fmt.Sprintf(", target_residues=%v", gm.targetResidues)
	}

	base += fmt.Sprintf(", mod_type='%s')", gm.globalModType)
	return base
}

// ToDict converts the modification to a dictionary representation
func (gm *GlobalModification) ToDict() map[string]interface{} {
	dict := gm.Modification.ToMap()
	dict["target_residues"] = gm.targetResidues
	dict["global_mod_type"] = gm.globalModType
	return dict
}
