package sequal

import (
	"crypto/sha256"
	"encoding/json"
	"fmt"
	"regexp"
	"strings"
)

// Modification represents a peptide modification with support for various ProForma 2.0 features
// including crosslinks, branches, ambiguity groups, and localization scores.
type Modification struct {
	BaseBlock
	source            *string
	originalValue     string
	crosslinkID       *string
	isCrosslinkRef    bool
	isBranchRef       bool
	isBranch          bool
	isAmbiguityRef    bool
	ambiguityGroup    *string
	inRange           bool
	rangeStart        *int
	rangeEnd          *int
	localizationScore *float64
	modValue          *ModificationValue

	regex       *regexp.Regexp
	modType     string
	labile      bool
	labilNumber int
	fullName    *string
	allFilled   bool
}

// KnownSources is a set of recognized modification source databases
var KnownSources = map[string]bool{
	"Unimod": true, "U": true, "PSI-MOD": true, "M": true,
	"RESID": true, "R": true, "XL-MOD": true, "X": true,
	"XLMOD": true, "GNO": true, "G": true, "MOD": true,
	"Obs": true, "Formula": true, "FORMULA": true, "GLYCAN": true,
	"Glycan": true, "Info": true, "INFO": true, "OBS": true,
	"XL": true,
}

// NewModification creates a new Modification instance with the specified parameters.
// It handles various ProForma 2.0 modification features including crosslinks, branches,
// ambiguity groups, and localization scores. The modType parameter must be one of the
// valid modification types (static, variable, terminal, ambiguous, crosslink, branch,
// gap, labile, unknown_position, global).
func NewModification(value string, position *int, regexPattern *string, fullName *string,
	modType string, labile bool, labilNumber int, mass float64, allFilled bool,
	crosslinkID *string, isCrosslinkRef bool, isBranchRef bool, isBranch bool,
	ambiguityGroup *string, isAmbiguityRef bool, inRange bool,
	rangeStart, rangeEnd *int, localizationScore *float64, modValue *ModificationValue) *Modification {

	var source *string
	originalValue := value

	if modValue == nil {
		valueWithCrosslink := value
		if crosslinkID != nil && !isCrosslinkRef {
			valueWithCrosslink = value + "#" + *crosslinkID
		}
		if ambiguityGroup != nil && !isAmbiguityRef {
			scoreStr := ""
			if localizationScore != nil {
				scoreStr = fmt.Sprintf("(%.2f)", *localizationScore)
			}
			valueWithCrosslink = value + "#" + *ambiguityGroup + scoreStr
		}
		modValue = NewModificationValue(valueWithCrosslink, &mass)
	}

	if len(value) > 0 && value[0] == '#' && isCrosslinkRef {
		clID := value[1:]
		crosslinkID = &clID
		value = "#" + clID
	}

	baseBlock := NewBaseBlock(value, position, true, &mass)

	validModTypes := map[string]bool{
		"static": true, "variable": true, "terminal": true, "ambiguous": true,
		"crosslink": true, "branch": true, "gap": true, "labile": true,
		"unknown_position": true, "global": true,
	}

	if (crosslinkID != nil || isCrosslinkRef) && modType != "crosslink" {
		modType = "crosslink"
	}

	if !validModTypes[modType] {
		validTypesList := make([]string, 0, len(validModTypes))
		for k := range validModTypes {
			validTypesList = append(validTypesList, k)
		}
		panic(fmt.Sprintf("mod_type must be one of: %s", strings.Join(validTypesList, ", ")))
	}

	var re *regexp.Regexp
	if regexPattern != nil {
		re = regexp.MustCompile(*regexPattern)
	}

	mod := &Modification{
		BaseBlock:         baseBlock,
		source:            source,
		originalValue:     originalValue,
		crosslinkID:       crosslinkID,
		isCrosslinkRef:    isCrosslinkRef,
		isBranchRef:       isBranchRef,
		isBranch:          isBranch,
		isAmbiguityRef:    isAmbiguityRef,
		ambiguityGroup:    ambiguityGroup,
		inRange:           inRange,
		rangeStart:        rangeStart,
		rangeEnd:          rangeEnd,
		localizationScore: localizationScore,
		modValue:          modValue,
		regex:             re,
		modType:           modType,
		labile:            labile,
		labilNumber:       labilNumber,
		fullName:          fullName,
		allFilled:         allFilled,
	}

	if modType == "labile" {
		mod.labile = true
	}
	if inRange {
		mod.modType = "ambiguous"
	}

	return mod
}

// GetValue returns the primary value of the modification.
// If a ModificationValue is set, it returns the primary value from that;
// otherwise, it returns the base block value.
func (m *Modification) GetValue() string {
	if m.modValue != nil {
		return m.modValue.GetPrimaryValue()
	}
	return m.BaseBlock.GetValue()
}

// GetMass returns the mass of the modification.
// If a ModificationValue is set, it returns the mass from that;
// otherwise, it returns the base block mass.
func (m *Modification) GetMass() *float64 {
	if m.modValue != nil {
		return m.modValue.GetMass()
	}
	return m.BaseBlock.GetMass()
}

// GetObservedMass returns the observed mass of the modification if available.
// This is only available through the ModificationValue.
func (m *Modification) GetObservedMass() *float64 {
	if m.modValue != nil {
		return m.modValue.GetObservedMass()
	}
	return nil
}

// GetAmbiguityGroup returns the ambiguity group identifier for ambiguous modifications.
func (m *Modification) GetAmbiguityGroup() *string {
	if m.modValue != nil {
		return m.modValue.GetAmbiguityGroup()
	}
	return nil
}

// IsAmbiguityRef returns true if this modification is a reference to an ambiguity group.
func (m *Modification) IsAmbiguityRef() bool {
	if m.modValue != nil {
		return m.modValue.IsAmbiguityRef()
	}
	return m.isAmbiguityRef
}

// GetSynonyms returns all synonyms of the modification value.
func (m *Modification) GetSynonyms() []string {
	return m.modValue.GetSynonyms()
}

// GetModificationValue returns the underlying ModificationValue object.
func (m *Modification) GetModificationValue() *ModificationValue {
	return m.modValue
}

// GetInfoTags returns the list of information tags associated with this modification.
func (m *Modification) GetInfoTags() []string {
	return m.modValue.GetInfoTags()
}

// GetCrosslinkID returns the crosslink identifier if this is a crosslink modification.
func (m *Modification) GetCrosslinkID() *string {
	if m.modValue != nil {
		return m.modValue.GetCrosslinkID()
	}
	return m.crosslinkID
}

// IsCrosslinkRef returns true if this modification is a reference to a crosslink.
func (m *Modification) IsCrosslinkRef() bool {
	if m.modValue != nil {
		return m.modValue.IsCrosslinkRef()
	}
	return m.isCrosslinkRef
}

// GetSource returns the modification database source (e.g., "Unimod", "PSI-MOD").
func (m *Modification) GetSource() *string {
	if m.modValue != nil {
		return m.modValue.GetSource()
	}
	return m.source
}

// GetOriginalValue returns the original modification value including any source prefix.
func (m *Modification) GetOriginalValue() string {
	return m.originalValue
}

// GetRegex returns the compiled regex pattern for finding modification sites in sequences.
func (m *Modification) GetRegex() *regexp.Regexp {
	return m.regex
}

// GetModType returns the modification type (e.g., "static", "variable", "terminal").
func (m *Modification) GetModType() string {
	return m.modType
}

// IsLabile returns true if the modification is labile (can be lost during fragmentation).
func (m *Modification) IsLabile() bool {
	return m.labile
}

// GetLabileNumber returns the labile fragmentation order number.
func (m *Modification) GetLabileNumber() int {
	return m.labilNumber
}

// GetFullName returns the full descriptive name of the modification if available.
func (m *Modification) GetFullName() *string {
	return m.fullName
}

// IsAllFilled returns true if the modification occurs at all expected sites.
func (m *Modification) IsAllFilled() bool {
	return m.allFilled
}

// FindPositions finds positions of the modification in the given sequence
func (m *Modification) FindPositions(seq string) [][]int {
	if m.regex == nil {
		panic(fmt.Sprintf("No regex pattern defined for modification '%s'", m.GetValue()))
	}

	var positions [][]int
	matches := m.regex.FindAllStringSubmatchIndex(seq, -1)

	for _, match := range matches {
		if len(match) > 2 { // Has groups
			for i := 0; i < len(match)/2; i++ {
				start := match[i*2]
				end := match[i*2+1]
				if start >= 0 && end >= 0 {
					positions = append(positions, []int{start, end})
				}
			}
		} else {
			positions = append(positions, []int{match[0], match[1]})
		}
	}

	return positions
}

// ToMap converts the modification to a map representation
func (m *Modification) ToMap() map[string]interface{} {
	result := m.BaseBlock.ToMap()

	var sourceStr *string
	if m.source != nil {
		sourceStr = m.source
	}

	var regexPattern *string
	if m.regex != nil {
		pattern := m.regex.String()
		regexPattern = &pattern
	}

	result["source"] = sourceStr
	result["original_value"] = m.originalValue
	result["regex_pattern"] = regexPattern
	result["full_name"] = m.fullName
	result["mod_type"] = m.modType
	result["labile"] = m.labile
	result["labile_number"] = m.labilNumber
	result["all_filled"] = m.allFilled
	result["crosslink_id"] = m.crosslinkID
	result["is_crosslink_ref"] = m.isCrosslinkRef

	return result
}

// Equal checks if two modifications are equal
func (m *Modification) Equal(other Modification) bool {
	if !m.BaseBlock.Equal(other.BaseBlock) {
		return false
	}

	mHash, err := m.Hash()
	if err != nil {
		return false
	}
	otherHash, err := other.Hash()
	if err != nil {
		return false
	}

	return mHash == otherHash
}

// Hash generates a hash for the modification
func (m *Modification) Hash() (string, error) {
	modMap := m.ToMap()
	jsonData, err := json.Marshal(modMap)
	if err != nil {
		return "", err
	}
	hash := sha256.Sum256(jsonData)
	return fmt.Sprintf("%x", hash), nil

}

// String returns a string representation of the modification
func (m *Modification) String() string {
	if m.isCrosslinkRef && m.crosslinkID != nil {
		return "#" + *m.crosslinkID
	}
	if m.isBranchRef {
		return "#BRANCH"
	}

	result := m.modValue.ToString()

	if m.crosslinkID != nil && !m.isCrosslinkRef {
		result += "#" + *m.crosslinkID
	}
	if m.isBranch && !m.isBranchRef {
		result += "#BRANCH"
	}
	if m.labile {
		result += fmt.Sprintf("%d", m.labilNumber)
	}

	return result
}

// HasAmbiguity checks if the modification has ambiguity
func (m *Modification) HasAmbiguity() bool {
	for _, v := range m.modValue.GetPipeValues() {
		if v.GetType() == PipeValueTypeAmbiguity {
			return true
		}
	}
	return false
}

// HasCrosslink checks if the modification has crosslink
func (m *Modification) HasCrosslink() bool {
	for _, v := range m.modValue.GetPipeValues() {
		if v.GetType() == PipeValueTypeCrosslink {
			return true
		}
	}
	return false
}

// HasBranch checks if the modification has branch
func (m *Modification) HasBranch() bool {
	for _, v := range m.modValue.GetPipeValues() {
		if v.GetType() == PipeValueTypeBranch {
			return true
		}
	}
	return false
}

// ToProforma converts the modification to ProForma notation string
func (m *Modification) ToProforma() string {
	if m.modValue != nil {
		seen := map[string]bool{}
		parts := []string{}

		for _, pv := range m.modValue.GetPipeValues() {
			modPart := ""
			if pv.GetSource() != nil {
				modPart = *pv.GetSource() + ":"
				if pv.GetMass() != nil {
					mass := *pv.GetMass()
					if mass > 0 {
						modPart += fmt.Sprintf("+%g", mass)
						seen[fmt.Sprintf("+%g", mass)] = true
					} else if mass < 0 {
						modPart += fmt.Sprintf("%g", mass)
						seen[fmt.Sprintf("%g", mass)] = true
					}
				} else {
					modPart += pv.GetValue()
				}
			} else {
				if pv.GetMass() != nil {
					mass := *pv.GetMass()
					if mass > 0 {
						modPart = fmt.Sprintf("+%g", mass)
					} else if mass < 0 {
						modPart = fmt.Sprintf("%g", mass)
					}
				} else if pv.GetType() == PipeValueTypeSynonym {
					modPart = pv.GetValue()
				} else {
					if !strings.Contains(pv.GetValue(), "#") {
						modPart = pv.GetValue()
					}
				}
			}

			if pv.GetType() == PipeValueTypeCrosslink && pv.GetCrosslinkID() != nil {
				modPart += "#" + *pv.GetCrosslinkID()
			} else if pv.GetType() == PipeValueTypeBranch && pv.IsBranch() {
				modPart += "#BRANCH"
			} else if pv.GetType() == PipeValueTypeAmbiguity && pv.GetAmbiguityGroup() != nil {
				scoreStr := ""
				if pv.GetLocalizationScore() != nil {
					scoreStr = fmt.Sprintf("(%.2f)", *pv.GetLocalizationScore())
				}
				modPart += "#" + *pv.GetAmbiguityGroup() + scoreStr
			}

			if _, exists := seen[modPart]; exists || modPart == "" {
				continue
			}
			parts = append(parts, modPart)
			seen[modPart] = true
		}

		return strings.Join(parts, "|")
	} else {
		if m.modValue != nil && m.BaseBlock.GetMass() != nil && (strings.HasPrefix(m.modValue.primaryValue, "+") || strings.HasPrefix(m.modValue.primaryValue, "-")) {
			return fmt.Sprintf("%g", *m.BaseBlock.GetMass())
		}
		if m.modValue != nil {
			return m.modValue.primaryValue
		}
		return m.BaseBlock.GetValue()
	}
}
