package sequal

import (
	"regexp"
	"strconv"
	"strings"
)

// PipeValueType represents the type of a pipe-separated value
type PipeValueType string

// Constants for different pipe value types
const (
	PipeValueTypeSynonym      PipeValueType = "synonym"
	PipeValueTypeInfoTag      PipeValueType = "info_tag"
	PipeValueTypeMass         PipeValueType = "mass"
	PipeValueTypeObservedMass PipeValueType = "observed_mass"
	PipeValueTypeCrosslink    PipeValueType = "crosslink"
	PipeValueTypeBranch       PipeValueType = "branch"
	PipeValueTypeAmbiguity    PipeValueType = "ambiguity"
	PipeValueTypeGlycan       PipeValueType = "glycan"
	PipeValueTypeGap          PipeValueType = "gap"
	PipeValueTypeFormula      PipeValueType = "formula"
)

// PipeValue represents a single pipe-separated value in a modification
type PipeValue struct {
	value             string
	valueType         PipeValueType
	crosslinkID       *string
	isBranch          bool
	isBranchRef       bool
	isCrosslinkRef    bool
	ambiguityGroup    *string
	isAmbiguityRef    bool
	localizationScore *float64
	source            *string
	originalValue     *string
	mass              *float64
	observedMass      *float64
	isValidGlycan     bool
	isValidFormula    bool
	assignedTypes     []PipeValueType
	charge            *string // ProForma 2.1: charge notation, e.g., "z+2", "z-1"
	chargeValue       *int    // ProForma 2.1: numeric charge value, e.g., 2, -1
}

// NewPipeValue creates a new PipeValue instance
func NewPipeValue(value string, valueType PipeValueType, originalValue string) *PipeValue {
	var origValue *string
	if originalValue != "" {
		origValue = &originalValue
	}

	pv := &PipeValue{
		value:          value,
		valueType:      valueType,
		originalValue:  origValue,
		isBranch:       false,
		isBranchRef:    false,
		isCrosslinkRef: false,
		isAmbiguityRef: false,
		isValidGlycan:  false,
		isValidFormula: false,
		assignedTypes:  make([]PipeValueType, 0),
	}

	pv.extractProperties()
	return pv
}

// extractProperties extracts special properties from the value based on type
func (pv *PipeValue) extractProperties() {
	// Handle crosslink values with #
	if pv.valueType == PipeValueTypeCrosslink && strings.Contains(pv.value, "#") {
		parts := strings.SplitN(pv.value, "#", 2)
		if parts[1] == "BRANCH" {
			pv.isBranch = true
		} else {
			crosslinkID := parts[1]
			pv.crosslinkID = &crosslinkID
		}
	} else if pv.valueType == PipeValueTypeAmbiguity && strings.Contains(pv.value, "#") {
		// Handle ambiguity values with #
		parts := strings.SplitN(pv.value, "#", 2)
		ambiguityGroup := parts[1]
		pv.ambiguityGroup = &ambiguityGroup

		// Check for localization score in parentheses
		if strings.Contains(ambiguityGroup, "(") && strings.Contains(ambiguityGroup, ")") {
			re := regexp.MustCompile(`\(([\d.]+)\)`)
			matches := re.FindStringSubmatch(ambiguityGroup)
			if len(matches) > 1 {
				if score, err := strconv.ParseFloat(matches[1], 64); err == nil {
					pv.localizationScore = &score
				}
			}
		}
	}

	// ProForma 2.1: Handle charged formulas (Section 11.1)
	// Format: Formula:value:z+2 or Formula:value:z-1
	if pv.valueType == PipeValueTypeFormula {
		// Check for charge notation pattern :z+N or :z-N
		re := regexp.MustCompile(`:z([+-]\d+)$`)
		matches := re.FindStringSubmatch(pv.value)
		if len(matches) > 1 {
			// Extract the charge string (e.g., "z+2")
			chargeStr := "z" + matches[1]
			pv.charge = &chargeStr

			// Extract the numeric charge value
			if chargeVal, err := strconv.Atoi(matches[1]); err == nil {
				pv.chargeValue = &chargeVal
			}

			// Remove the charge notation from the value
			pv.value = re.ReplaceAllString(pv.value, "")
		}
	}
}

// String returns the string representation
func (pv *PipeValue) String() string {
	return pv.value
}

// GetType returns the value type
func (pv *PipeValue) GetType() PipeValueType {
	return pv.valueType
}

// SetType changes the value type
func (pv *PipeValue) SetType(valueType PipeValueType) {
	pv.valueType = valueType
	if len(pv.assignedTypes) > 0 {
		pv.assignedTypes[0] = valueType
	} else {
		pv.AssignType(valueType)
	}
}

// AssignType adds a type to the assigned types list
func (pv *PipeValue) AssignType(valueType PipeValueType) {
	// Check if this type is already assigned
	for _, t := range pv.assignedTypes {
		if t == valueType {
			return
		}
	}
	pv.assignedTypes = append(pv.assignedTypes, valueType)
}

// GetValue returns the value
func (pv *PipeValue) GetValue() string {
	return pv.value
}

// GetCrosslinkID returns the crosslink ID
func (pv *PipeValue) GetCrosslinkID() *string {
	return pv.crosslinkID
}

// IsBranch returns whether this is a branch
func (pv *PipeValue) IsBranch() bool {
	return pv.isBranch
}

// IsBranchRef returns whether this is a branch reference
func (pv *PipeValue) IsBranchRef() bool {
	return pv.isBranchRef
}

// IsCrosslinkRef returns whether this is a crosslink reference
func (pv *PipeValue) IsCrosslinkRef() bool {
	return pv.isCrosslinkRef
}

// GetAmbiguityGroup returns the ambiguity group
func (pv *PipeValue) GetAmbiguityGroup() *string {
	return pv.ambiguityGroup
}

// IsAmbiguityRef returns whether this is an ambiguity reference
func (pv *PipeValue) IsAmbiguityRef() bool {
	return pv.isAmbiguityRef
}

// GetLocalizationScore returns the localization score
func (pv *PipeValue) GetLocalizationScore() *float64 {
	return pv.localizationScore
}

// GetSource returns the source
func (pv *PipeValue) GetSource() *string {
	return pv.source
}

// SetSource sets the source
func (pv *PipeValue) SetSource(source string) {
	pv.source = &source
}

// GetMass returns the mass
func (pv *PipeValue) GetMass() *float64 {
	return pv.mass
}

// SetMass sets the mass
func (pv *PipeValue) SetMass(mass float64) {
	pv.mass = &mass
}

// GetObservedMass returns the observed mass
func (pv *PipeValue) GetObservedMass() *float64 {
	return pv.observedMass
}

// SetObservedMass sets the observed mass
func (pv *PipeValue) SetObservedMass(mass float64) {
	pv.observedMass = &mass
}

// GetCharge returns the charge notation (ProForma 2.1)
func (pv *PipeValue) GetCharge() *string {
	return pv.charge
}

// GetChargeValue returns the numeric charge value (ProForma 2.1)
func (pv *PipeValue) GetChargeValue() *int {
	return pv.chargeValue
}

// IsValidGlycan returns whether this is a valid glycan
func (pv *PipeValue) IsValidGlycan() bool {
	return pv.isValidGlycan
}

// IsValidFormula returns whether this is a valid formula
func (pv *PipeValue) IsValidFormula() bool {
	return pv.isValidFormula
}
