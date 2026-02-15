package sequal

import (
	"regexp"
	"sort"
	"strconv"
	"strings"
	"unicode"
)

// ModificationValue represents a modification value with unified pipe value handling
type ModificationValue struct {
	primaryValue string
	source       *string
	mass         *float64
	pipeValues   []*PipeValue
	knownSources map[string]bool
}

// NewModificationValue creates a new ModificationValue instance
func NewModificationValue(value string, mass *float64) *ModificationValue {
	mv := &ModificationValue{
		pipeValues: make([]*PipeValue, 0),
		knownSources: map[string]bool{
			"Unimod": true, "U": true, "PSI-MOD": true, "M": true,
			"RESID": true, "R": true, "XL-MOD": true, "X": true,
			"XLMOD": true, "GNO": true, "G": true, "MOD": true,
			"Obs": true, "Formula": true, "FORMULA": true, "GLYCAN": true,
			"Glycan": true, "Info": true, "OBS": true,
			"INFO": true, "XL": true,
		},
		mass: mass,
	}
	mv._parseValue(value)
	return mv
}

// _parseValue parses modification value with unified pipe value handling
func (mv *ModificationValue) _parseValue(value string) {
	if strings.Contains(value, "|") {
		components := strings.Split(value, "|")
		mv.processPrimaryValue(components[0])
		for _, component := range components[1:] {
			mv._processPipeComponent(component)
		}
	} else {
		mv.processPrimaryValue(value)
	}
}

// GetSource returns the source of the modification value
func (mv *ModificationValue) GetSource() *string {
	return mv.source
}

// GetPrimaryValue returns the primary value
func (mv *ModificationValue) GetPrimaryValue() string {
	return mv.primaryValue
}

// GetMass returns the mass value
func (mv *ModificationValue) GetMass() *float64 {
	return mv.mass
}

// GetSynonyms returns list of synonyms
func (mv *ModificationValue) GetSynonyms() []string {
	result := make([]string, 0)
	for _, pv := range mv.pipeValues {
		if pv.valueType == PipeValueTypeSynonym {
			result = append(result, pv.value)
		}
	}
	return result
}

// GetObservedMass returns the observed mass if any
func (mv *ModificationValue) GetObservedMass() *float64 {
	for _, pv := range mv.pipeValues {
		if pv.valueType == PipeValueTypeObservedMass && pv.observedMass != nil {
			return pv.observedMass
		}
	}
	return nil
}

// GetPipeValues returns all pipe values
func (mv *ModificationValue) GetPipeValues() []*PipeValue {
	return mv.pipeValues
}

// GetInfoTags returns list of info tags
func (mv *ModificationValue) GetInfoTags() []string {
	result := make([]string, 0)
	for _, pv := range mv.pipeValues {
		if pv.valueType == PipeValueTypeInfoTag {
			result = append(result, pv.value)
		}
	}
	return result
}

// GetCrosslinkID returns the crosslink ID if any
func (mv *ModificationValue) GetCrosslinkID() *string {
	for _, pv := range mv.pipeValues {
		if pv.valueType == PipeValueTypeCrosslink && pv.crosslinkID != nil {
			return pv.crosslinkID
		}
	}
	return nil
}

// IsBranch checks if modification is a branch
func (mv *ModificationValue) IsBranch() bool {
	for _, pv := range mv.pipeValues {
		if pv.valueType == PipeValueTypeBranch && pv.isBranch {
			return true
		}
	}
	return false
}

// IsBranchRef checks if modification is a branch reference
func (mv *ModificationValue) IsBranchRef() bool {
	for _, pv := range mv.pipeValues {
		if pv.valueType == PipeValueTypeBranch && pv.isBranchRef {
			return true
		}
	}
	return false
}

// IsCrosslinkRef checks if modification is a crosslink reference
func (mv *ModificationValue) IsCrosslinkRef() bool {
	for _, pv := range mv.pipeValues {
		if pv.valueType == PipeValueTypeCrosslink && pv.isCrosslinkRef {
			return true
		}
	}
	return false
}

// GetAmbiguityGroup returns the ambiguity group if any
func (mv *ModificationValue) GetAmbiguityGroup() *string {
	for _, pv := range mv.pipeValues {
		if pv.valueType == PipeValueTypeAmbiguity && pv.ambiguityGroup != nil {
			return pv.ambiguityGroup
		}
	}
	return nil
}

// IsAmbiguityRef checks if modification is an ambiguity reference
func (mv *ModificationValue) IsAmbiguityRef() bool {
	for _, pv := range mv.pipeValues {
		if pv.valueType == PipeValueTypeAmbiguity && pv.isAmbiguityRef {
			return true
		}
	}
	return false
}

// ToString returns string representation
func (mv *ModificationValue) ToString() string {
	parts := make([]string, 0)
	seen := make(map[string]bool)

	for _, pv := range mv.pipeValues {
		if seen[pv.value] {
			continue
		}
		parts = append(parts, pv.value)
		seen[pv.value] = true
	}

	return strings.Join(parts, "|")
}

// _validateFormula statically validates a chemical formula
func _validateFormula(formula string) bool {
	// Empty formula is invalid
	if len(strings.TrimSpace(formula)) == 0 {
		return false
	}

	if strings.Count(formula, "[") != strings.Count(formula, "]") {
		return false
	}

	formulaNoSpaces := strings.ReplaceAll(formula, " ", "")

	i := 0
	for i < len(formulaNoSpaces) {
		if formulaNoSpaces[i] == '[' {
			endBracket := strings.Index(formulaNoSpaces[i:], "]")
			if endBracket == -1 {
				return false
			}
			endBracket += i

			isotopePart := formulaNoSpaces[i+1 : endBracket]
			matched, _ := regexp.MatchString(`\d+[A-Z][a-z]?(-?\d+)?`, isotopePart)
			if !matched {
				return false
			}
			i = endBracket + 1

			if i < len(formulaNoSpaces) && (formulaNoSpaces[i] == '-' || unicode.IsDigit(rune(formulaNoSpaces[i]))) {
				j := i
				if formulaNoSpaces[j] == '-' {
					j++
				}
				for j < len(formulaNoSpaces) && unicode.IsDigit(rune(formulaNoSpaces[j])) {
					j++
				}
				i = j
			}
		} else if unicode.IsUpper(rune(formulaNoSpaces[i])) {
			i++
			if i < len(formulaNoSpaces) && unicode.IsLower(rune(formulaNoSpaces[i])) {
				i++
			}

			if i < len(formulaNoSpaces) && (formulaNoSpaces[i] == '-' || unicode.IsDigit(rune(formulaNoSpaces[i]))) {
				j := i
				if formulaNoSpaces[j] == '-' {
					j++
				}
				for j < len(formulaNoSpaces) && unicode.IsDigit(rune(formulaNoSpaces[j])) {
					j++
				}
				i = j
			}
		} else {
			return false
		}
	}

	return true
}

// _validateGlycan statically validates a glycan string
func _validateGlycan(glycan string) bool {
	glycanClean := strings.ReplaceAll(glycan, " ", "")

	monos := []string{
		"Hex", "HexNAc", "HexS", "HexP", "HexNAcS",
		"dHex", "NeuAc", "NeuGc", "Pen", "Fuc",
	}

	// Sort by length (longest first) to avoid partial matches
	sort.Slice(monos, func(i, j int) bool {
		return len(monos[i]) > len(monos[j])
	})

	monoPattern := "^("
	for i, mono := range monos {
		if i > 0 {
			monoPattern += "|"
		}
		monoPattern += regexp.QuoteMeta(mono)
	}
	monoPattern += `)((\\([1-9]\\d*\\))|[1-9]\\d*)?`

	re := regexp.MustCompile(monoPattern)

	i := 0
	for i < len(glycanClean) {
		if glycanClean[i] == '{' {
			closeBrace := strings.Index(glycanClean[i:], "}")
			if closeBrace == -1 {
				return false
			}
			closeBrace += i

			i = closeBrace + 1
			isAtEnd := i == len(glycanClean)

			if i < len(glycanClean) && glycanClean[i] == '(' {
				closeParen := strings.Index(glycanClean[i:], ")")
				if closeParen == -1 {
					return false
				}
				closeParen += i
				countStr := glycanClean[i+1 : closeParen]
				matched, _ := regexp.MatchString(`^[1-9]\d*$`, countStr)
				if !matched {
					return false
				}
				i = closeParen + 1
			} else if i < len(glycanClean) && unicode.IsDigit(rune(glycanClean[i])) {
				start := i
				for i < len(glycanClean) && unicode.IsDigit(rune(glycanClean[i])) {
					i++
				}
				countStr := glycanClean[start:i]
				matched, _ := regexp.MatchString(`^[1-9]\d*$`, countStr)
				if !matched {
					return false
				}
			} else if !isAtEnd {
				return false
			}
			continue
		}

		match := re.FindStringSubmatch(glycanClean[i:])
		if match == nil {
			return false
		}

		monoLength := len(match[0])
		hasCount := len(match) > 2 && match[2] != ""
		i += monoLength

		isAtEnd := i == len(glycanClean)
		if !hasCount && !isAtEnd {
			return false
		}
	}

	return i == len(glycanClean)
}

// Get returns pipe value at index
func (mv *ModificationValue) Get(index int) *PipeValue {
	if index >= 0 && index < len(mv.pipeValues) {
		return mv.pipeValues[index]
	}
	return nil
}

// Len returns number of pipe values
func (mv *ModificationValue) Len() int {
	return len(mv.pipeValues)
}

// ForEach iterates through all pipe values
func (mv *ModificationValue) ForEach(fn func(pv *PipeValue)) {
	for _, pv := range mv.pipeValues {
		fn(pv)
	}
}

// processPrimaryValue processes the primary value component
func (mv *ModificationValue) processPrimaryValue(value string) {
	// Handle branch reference
	if value == "#BRANCH" {
		mv.primaryValue = ""
		pipeVal := NewPipeValue(value, PipeValueTypeBranch, value)
		pipeVal.isBranchRef = true
		pipeVal.isBranch = true
		mv.pipeValues = append(mv.pipeValues, pipeVal)
		return
	} else if strings.HasPrefix(value, "#") {
		// Handle crosslink or ambiguity reference
		mv.primaryValue = ""
		valueType := PipeValueTypeCrosslink
		if strings.Contains(value[1:], "(") && strings.Contains(value[1:], ")") {
			valueType = PipeValueTypeAmbiguity
		}

		pipeVal := NewPipeValue(value, valueType, value)
		pipeVal.isCrosslinkRef = valueType == PipeValueTypeCrosslink
		pipeVal.isAmbiguityRef = valueType == PipeValueTypeAmbiguity

		if pipeVal.isCrosslinkRef {
			crosslinkID := value[1:]
			pipeVal.crosslinkID = &crosslinkID
		}

		if pipeVal.isAmbiguityRef {
			ambiguityGroup := value[1:]
			pipeVal.ambiguityGroup = &ambiguityGroup

			if strings.Contains(ambiguityGroup, "(") && strings.Contains(ambiguityGroup, ")") {
				re := regexp.MustCompile(`\(([\d.]+)\)`)
				matches := re.FindStringSubmatch(ambiguityGroup)
				if len(matches) > 1 {
					score, err := strconv.ParseFloat(matches[1], 64)
					if err == nil {
						pipeVal.localizationScore = &score
						// Remove the score part from ambiguity group
						cleanGroup := re.ReplaceAllString(ambiguityGroup, "")
						pipeVal.ambiguityGroup = &cleanGroup
					}
				}
			}
		}

		mv.pipeValues = append(mv.pipeValues, pipeVal)
		return
	}

	// Handle source prefix
	if strings.Contains(value, ":") {
		parts := strings.SplitN(value, ":", 2)
		if mv.knownSources[parts[0]] {
			source := parts[0]
			mv.source = &source

			// Handle empty value after source
			if parts[1] == "" {
				mv.primaryValue = ""
				return
			}

			valueStr := parts[1]
			mv.primaryValue = valueStr

			// ProForma 2.1: Extract charge notation for formulas
			var chargeStr *string
			var chargeValue *int
			if strings.ToUpper(source) == "FORMULA" {
				// Check for charge notation pattern :z+N or :z-N
				re := regexp.MustCompile(`:z([+-]\d+)$`)
				matches := re.FindStringSubmatch(valueStr)
				if len(matches) > 1 {
					// Extract the charge string (e.g., "z+2")
					charge := "z" + matches[1]
					chargeStr = &charge

					// Extract the numeric charge value
					if cVal, err := strconv.Atoi(matches[1]); err == nil {
						chargeValue = &cVal
					}

					// Remove the charge notation from the value
					valueStr = re.ReplaceAllString(valueStr, "")
					mv.primaryValue = valueStr
				}
			}

			// Handle special cases with # symbol
			if strings.Contains(valueStr, "#") {
				valueParts := strings.SplitN(valueStr, "#", 2)
				baseValue := valueParts[0]
				specialPart := valueParts[1]

				// Create pipe value for base value
				pipeVal := NewPipeValue(baseValue, PipeValueTypeSynonym, valueStr)
				pipeVal.source = &source
				mv.pipeValues = append(mv.pipeValues, pipeVal)

				// Handle branch, ambiguity, or crosslink
				if specialPart == "BRANCH" {
					branchVal := NewPipeValue(valueStr, PipeValueTypeBranch, valueStr)
					branchVal.isBranch = true
					branchVal.source = &source
					mv.pipeValues = append(mv.pipeValues, branchVal)
				} else if strings.Contains(specialPart, "(") && strings.Contains(specialPart, ")") {
					ambVal := NewPipeValue(valueStr, PipeValueTypeAmbiguity, valueStr)
					ambiguityGroup := specialPart
					ambVal.ambiguityGroup = &ambiguityGroup
					ambVal.source = &source

					re := regexp.MustCompile(`\(([\d.]+)\)`)
					matches := re.FindStringSubmatch(specialPart)
					if len(matches) > 1 {
						score, err := strconv.ParseFloat(matches[1], 64)
						if err == nil {
							ambVal.localizationScore = &score
							// Remove the score part from ambiguity group
							cleanGroup := re.ReplaceAllString(ambiguityGroup, "")
							ambVal.ambiguityGroup = &cleanGroup
						}
					}

					mv.pipeValues = append(mv.pipeValues, ambVal)
				} else {
					xlVal := NewPipeValue(valueStr, PipeValueTypeCrosslink, valueStr)
					xlVal.crosslinkID = &specialPart
					xlVal.source = &source
					mv.pipeValues = append(mv.pipeValues, xlVal)
				}
			} else {
				// Simple value with source
				pipeVal := NewPipeValue(valueStr, PipeValueTypeSynonym, valueStr)
				pipeVal.source = &source
				// ProForma 2.1: Set charge if present (for formulas)
				if strings.ToUpper(source) == "FORMULA" {
					pipeVal.SetType(PipeValueTypeFormula)
					pipeVal.isValidFormula = validateFormula(valueStr)
					pipeVal.charge = chargeStr
					pipeVal.chargeValue = chargeValue
				} else if strings.ToUpper(source) == "GLYCAN" {
					// ProForma 2.1: Validate glycan (including custom monosaccharides)
					pipeVal.SetType(PipeValueTypeGlycan)
					pipeVal.isValidGlycan = validateGlycan(valueStr)
				}
				mv.pipeValues = append(mv.pipeValues, pipeVal)
			}
		} else if strings.ToUpper(parts[0]) == "MASS" {
			// Mass specification
			massStr := parts[1]
			mv.primaryValue = massStr

			mass, err := strconv.ParseFloat(massStr, 64)
			if err == nil {
				mv.mass = &mass
				massVal := NewPipeValue(massStr, PipeValueTypeMass, value)
				massVal.mass = &mass
				mv.pipeValues = append(mv.pipeValues, massVal)
			} else {
				infoVal := NewPipeValue(value, PipeValueTypeInfoTag, value)
				mv.pipeValues = append(mv.pipeValues, infoVal)
			}
		} else {
			// Unknown source, treat as info tag
			mv.primaryValue = value
			infoVal := NewPipeValue(value, PipeValueTypeInfoTag, value)
			mv.pipeValues = append(mv.pipeValues, infoVal)
		}
	} else {
		// No source prefix
		if strings.Contains(value, "#") {
			valueParts := strings.SplitN(value, "#", 2)
			baseValue := valueParts[0]
			specialPart := valueParts[1]

			mv.primaryValue = baseValue

			// Base value as synonym
			if baseValue != "" {
				synVal := NewPipeValue(baseValue, PipeValueTypeSynonym, value)
				mv.pipeValues = append(mv.pipeValues, synVal)
			}

			// Handle special part
			if specialPart == "BRANCH" {
				branchVal := NewPipeValue(value, PipeValueTypeBranch, value)
				branchVal.isBranch = true
				mv.pipeValues = append(mv.pipeValues, branchVal)
			} else if strings.Contains(specialPart, "(") && strings.Contains(specialPart, ")") {
				ambVal := NewPipeValue(value, PipeValueTypeAmbiguity, value)
				ambiguityGroup := specialPart
				ambVal.ambiguityGroup = &ambiguityGroup

				re := regexp.MustCompile(`\(([\d.]+)\)`)
				matches := re.FindStringSubmatch(specialPart)
				if len(matches) > 1 {
					score, err := strconv.ParseFloat(matches[1], 64)
					if err == nil {
						ambVal.localizationScore = &score
						// Remove the score part from ambiguity group
						cleanGroup := re.ReplaceAllString(ambiguityGroup, "")
						ambVal.ambiguityGroup = &cleanGroup
					}
				}

				mv.pipeValues = append(mv.pipeValues, ambVal)
			} else {
				xlVal := NewPipeValue(value, PipeValueTypeCrosslink, value)
				xlVal.crosslinkID = &specialPart
				mv.pipeValues = append(mv.pipeValues, xlVal)
			}
		} else {
			// Regular value (could be mass or synonym)
			mv.primaryValue = value

			// Try to parse as mass
			if strings.HasPrefix(value, "+") || strings.HasPrefix(value, "-") {
				massStr := value
				if strings.HasPrefix(value, "+") {
					massStr = value[1:] // Remove + sign
				}

				mass, err := strconv.ParseFloat(massStr, 64)
				if err == nil {
					if strings.HasPrefix(value, "-") {
						mass = -mass
					}

					mv.mass = &mass
					massVal := NewPipeValue(value, PipeValueTypeMass, value)
					massVal.mass = &mass
					mv.pipeValues = append(mv.pipeValues, massVal)
					return
				}
			}

			// Not a mass, treat as synonym
			synVal := NewPipeValue(value, PipeValueTypeSynonym, value)
			mv.pipeValues = append(mv.pipeValues, synVal)

			// Try to validate as glycan or formula
			if validateGlycan(value) {
				glycanVal := NewPipeValue(value, PipeValueTypeGlycan, value)
				glycanVal.isValidGlycan = true
				mv.pipeValues = append(mv.pipeValues, glycanVal)
			} else if validateFormula(value) {
				formulaVal := NewPipeValue(value, PipeValueTypeFormula, value)
				formulaVal.isValidFormula = true
				mv.pipeValues = append(mv.pipeValues, formulaVal)
			}
		}
	}
}

// processPipeComponent processes a single pipe-separated component
func (mv *ModificationValue) _processPipeComponent(component string) {
	if component == "#BRANCH" {
		pipeVal := NewPipeValue(component, PipeValueTypeBranch, component)
		pipeVal.isBranchRef = true
		pipeVal.isBranch = true
		mv.pipeValues = append(mv.pipeValues, pipeVal)
		return
	} else if strings.HasPrefix(component, "#") {
		var pipeValType PipeValueType
		if strings.HasPrefix(component[1:], "XL") {
			pipeValType = PipeValueTypeCrosslink
		} else {
			pipeValType = PipeValueTypeAmbiguity
		}

		pipeVal := NewPipeValue(component, pipeValType, component)
		pipeVal.isCrosslinkRef = pipeVal.valueType == PipeValueTypeCrosslink
		pipeVal.isAmbiguityRef = pipeVal.valueType == PipeValueTypeAmbiguity

		if pipeVal.isCrosslinkRef {
			crosslinkID := component[1:]
			pipeVal.crosslinkID = &crosslinkID
		}

		if pipeVal.isAmbiguityRef {
			ambiguityGroup := component[1:]
			pipeVal.ambiguityGroup = &ambiguityGroup

			if strings.Contains(ambiguityGroup, "(") && strings.Contains(ambiguityGroup, ")") {
				re := regexp.MustCompile(`\(([\d.]+)\)`)
				matches := re.FindStringSubmatch(ambiguityGroup)
				if len(matches) > 1 {
					if score, err := strconv.ParseFloat(matches[1], 64); err == nil {
						pipeVal.localizationScore = &score
						// Remove the score part from ambiguity group
						cleanGroup := re.ReplaceAllString(ambiguityGroup, "")
						pipeVal.ambiguityGroup = &cleanGroup
					}
				}
			}
		}

		mv.pipeValues = append(mv.pipeValues, pipeVal)
		return
	}

	if strings.Contains(component, ":") {
		parts := strings.SplitN(component, ":", 2)
		if isInKnownSources(parts[0], mv.knownSources) {
			source := parts[0]
			value := parts[1]

			// ProForma 2.1: Extract charge notation for formulas
			var chargeStr *string
			var chargeValue *int
			if strings.ToUpper(source) == "FORMULA" {
				// Check for charge notation pattern :z+N or :z-N
				re := regexp.MustCompile(`:z([+-]\d+)$`)
				matches := re.FindStringSubmatch(value)
				if len(matches) > 1 {
					// Extract the charge string (e.g., "z+2")
					charge := "z" + matches[1]
					chargeStr = &charge

					// Extract the numeric charge value
					if cVal, err := strconv.Atoi(matches[1]); err == nil {
						chargeValue = &cVal
					}

					// Remove the charge notation from the value
					value = re.ReplaceAllString(value, "")
				}
			}

			if strings.Contains(value, "#") {
				pvParts := strings.SplitN(value, "#", 2)
				value = pvParts[0]
				isValidGlycan := false
				isValidFormula := false

				if strings.ToUpper(source) == "FORMULA" {
					isValidFormula = validateFormula(value)
				} else if strings.ToUpper(source) == "GLYCAN" {
					isValidGlycan = validateGlycan(value)
				}

				var pipeVal *PipeValue

				if source == "XL" || source == "XLMOD" || source == "XL-MOD" || source == "X" {
					pipeVal = NewPipeValue(value, PipeValueTypeCrosslink, component)
					pipeVal.source = &source
					crosslinkID := pvParts[1]
					pipeVal.crosslinkID = &crosslinkID
				} else if pvParts[1] == "BRANCH" {
					pipeVal = NewPipeValue(value, PipeValueTypeBranch, component)
					pipeVal.source = &source
					pipeVal.isBranch = true
				} else if strings.ToUpper(source) == "GLYCAN" {
					pipeVal = NewPipeValue(value, PipeValueTypeGlycan, component)
					pipeVal.source = &source
					pipeVal.isValidGlycan = isValidGlycan
				} else if strings.ToUpper(source) == "GNO" || strings.ToUpper(*mv.source) == "G" {
					pipeVal = NewPipeValue(value, PipeValueTypeGlycan, component)
					pipeVal.source = &source
					pipeVal.isValidGlycan = true
				} else if strings.ToUpper(source) == "FORMULA" {
					pipeVal = NewPipeValue(value, PipeValueTypeFormula, component)
					pipeVal.source = &source
					pipeVal.isValidFormula = isValidFormula
					// ProForma 2.1: Set charge if present
					pipeVal.charge = chargeStr
					pipeVal.chargeValue = chargeValue
				} else {
					pipeVal = NewPipeValue(value, PipeValueTypeAmbiguity, component)
					pipeVal.source = &source

					if isValidGlycan {
						pipeVal.AssignType(PipeValueTypeGlycan)
						pipeVal.isValidGlycan = true
					} else if isValidFormula {
						pipeVal.AssignType(PipeValueTypeFormula)
						pipeVal.isValidFormula = true
					}

					if strings.ToUpper(source) == "GNO" || strings.ToUpper(*mv.source) == "G" {
						pipeVal.AssignType(PipeValueTypeGap)
					}

					ambiguityGroup := pvParts[1]
					pipeVal.ambiguityGroup = &ambiguityGroup

					if strings.Contains(ambiguityGroup, "(") && strings.Contains(ambiguityGroup, ")") {
						re := regexp.MustCompile(`\(([\d.]+)\)`)
						matches := re.FindStringSubmatch(ambiguityGroup)
						if len(matches) > 1 {
							if score, err := strconv.ParseFloat(matches[1], 64); err == nil {
								pipeVal.localizationScore = &score
								// Remove the score part from ambiguity group
								cleanGroup := re.ReplaceAllString(ambiguityGroup, "")
								pipeVal.ambiguityGroup = &cleanGroup
							}
						}
					}
				}

				mv.pipeValues = append(mv.pipeValues, pipeVal)

			} else {
				var pipeVal *PipeValue

				if strings.ToUpper(source) == "INFO" {
					pipeVal = NewPipeValue(value, PipeValueTypeInfoTag, component)
				} else if strings.ToUpper(source) == "OBS" {
					pipeVal = NewPipeValue(value, PipeValueTypeObservedMass, component)
					if observedMass, err := strconv.ParseFloat(value, 64); err == nil {
						pipeVal.observedMass = &observedMass
					}
				} else if strings.ToUpper(source) == "GLYCAN" {
					isValidGlycan := validateGlycan(value)
					pipeVal = NewPipeValue(value, PipeValueTypeGlycan, component)
					pipeVal.isValidGlycan = isValidGlycan
				} else if strings.ToUpper(source) == "GNO" || strings.ToUpper(source) == "G" {
					pipeVal = NewPipeValue(value, PipeValueTypeGlycan, component)
					pipeVal.isValidGlycan = true
				} else if strings.ToUpper(source) == "FORMULA" {
					isValidFormula := validateFormula(value)
					pipeVal = NewPipeValue(value, PipeValueTypeFormula, component)
					pipeVal.isValidFormula = isValidFormula
					// ProForma 2.1: Set charge if present
					pipeVal.charge = chargeStr
					pipeVal.chargeValue = chargeValue
				} else {
					pipeVal = NewPipeValue(value, PipeValueTypeSynonym, component)
				}

				pipeVal.source = &source
				mv.pipeValues = append(mv.pipeValues, pipeVal)
			}

		} else if strings.ToUpper(parts[0]) == "MASS" {
			if mass, err := strconv.ParseFloat(parts[1], 64); err == nil {
				pipeVal := NewPipeValue(parts[1], PipeValueTypeMass, component)
				pipeVal.mass = &mass

				if strings.Contains(parts[1], "#") {
					pvParts := strings.SplitN(parts[1], "#", 2)
					pipeVal.value = pvParts[0]

					if pvParts[1] == "BRANCH" {
						pipeVal.isBranch = true
						pipeVal.AssignType(PipeValueTypeBranch)
					} else if strings.HasPrefix(pvParts[1], "XL") {
						crosslinkID := pvParts[1]
						pipeVal.crosslinkID = &crosslinkID
						pipeVal.AssignType(PipeValueTypeCrosslink)
					} else {
						ambiguityGroup := pvParts[1]
						pipeVal.ambiguityGroup = &ambiguityGroup
						pipeVal.AssignType(PipeValueTypeAmbiguity)

						if strings.Contains(ambiguityGroup, "(") && strings.Contains(ambiguityGroup, ")") {
							re := regexp.MustCompile(`\(([\d.]+)\)`)
							matches := re.FindStringSubmatch(ambiguityGroup)
							if len(matches) > 1 {
								if score, err := strconv.ParseFloat(matches[1], 64); err == nil {
									pipeVal.localizationScore = &score
									// Remove the score part from ambiguity group
									cleanGroup := re.ReplaceAllString(ambiguityGroup, "")
									pipeVal.ambiguityGroup = &cleanGroup
								}
							}
						}
					}

					pipeVal.AssignType(PipeValueTypeMass)
				}

				mv.pipeValues = append(mv.pipeValues, pipeVal)

			} else {
				mv.pipeValues = append(mv.pipeValues, NewPipeValue(component, PipeValueTypeSynonym, component))
			}

		} else {
			mv.pipeValues = append(mv.pipeValues, NewPipeValue(component, PipeValueTypeSynonym, component))
		}

	} else {
		if strings.Contains(component, "#") {
			parts := strings.SplitN(component, "#", 2)
			value := parts[0]

			var pipeVal *PipeValue

			if parts[1] == "BRANCH" {
				pipeVal = NewPipeValue(value, PipeValueTypeBranch, component)
				pipeVal.isBranch = true
			} else if strings.HasPrefix(parts[1], "XL") {
				pipeVal = NewPipeValue(value, PipeValueTypeCrosslink, component)
				crosslinkID := parts[1]
				pipeVal.crosslinkID = &crosslinkID
			} else {
				pipeVal = NewPipeValue(value, PipeValueTypeAmbiguity, component)
				ambiguityGroup := parts[1]
				pipeVal.ambiguityGroup = &ambiguityGroup

				if strings.Contains(ambiguityGroup, "(") && strings.Contains(ambiguityGroup, ")") {
					re := regexp.MustCompile(`\(([\d.]+)\)`)
					matches := re.FindStringSubmatch(ambiguityGroup)
					if len(matches) > 1 {
						if score, err := strconv.ParseFloat(matches[1], 64); err == nil {
							pipeVal.localizationScore = &score
							// Remove the score part from ambiguity group
							cleanGroup := re.ReplaceAllString(ambiguityGroup, "")
							pipeVal.ambiguityGroup = &cleanGroup
						}
					}
				}
			}

			if (strings.HasPrefix(value, "+") || strings.HasPrefix(value, "-")) && containsDigit(value) {
				if mass, err := strconv.ParseFloat(value, 64); err == nil {
					pipeVal.mass = &mass
					pipeVal.AssignType(PipeValueTypeMass)
				}
			} else {
				pipeVal.AssignType(PipeValueTypeSynonym)
			}

			mv.pipeValues = append(mv.pipeValues, pipeVal)

		} else {
			if (strings.HasPrefix(component, "+") || strings.HasPrefix(component, "-")) && containsDigit(component) {
				if mass, err := strconv.ParseFloat(component, 64); err == nil {
					pipeVal := NewPipeValue(component, PipeValueTypeMass, component)
					pipeVal.mass = &mass
					mv.pipeValues = append(mv.pipeValues, pipeVal)
				} else {
					mv.pipeValues = append(mv.pipeValues, NewPipeValue(component, PipeValueTypeSynonym, component))
				}
			} else {
				mv.pipeValues = append(mv.pipeValues, NewPipeValue(component, PipeValueTypeSynonym, component))
			}
		}
	}
}

func containsDigit(s string) bool {
	for _, c := range s {
		if unicode.IsDigit(c) {
			return true
		}
	}
	return false
}

func isInKnownSources(source string, knownSources map[string]bool) bool {
	_, ok := knownSources[source]
	return ok
}

// validateFormula statically validates a chemical formula
func validateFormula(formula string) bool {
	// Empty formula is invalid
	if strings.TrimSpace(formula) == "" {
		return false
	}

	if strings.Count(formula, "[") != strings.Count(formula, "]") {
		return false
	}

	formulaNoSpaces := strings.ReplaceAll(formula, " ", "")

	i := 0
	for i < len(formulaNoSpaces) {
		if formulaNoSpaces[i] == '[' {
			endBracket := strings.Index(formulaNoSpaces[i:], "]")
			if endBracket == -1 {
				return false
			}
			endBracket += i

			isotopePart := formulaNoSpaces[i+1 : endBracket]
			matched, _ := regexp.MatchString(`\d+[A-Z][a-z]?(-?\d+)?`, isotopePart)
			if !matched {
				return false
			}
			i = endBracket + 1

			if i < len(formulaNoSpaces) && (formulaNoSpaces[i] == '-' || unicode.IsDigit(rune(formulaNoSpaces[i]))) {
				j := i
				if formulaNoSpaces[j] == '-' {
					j++
				}
				for j < len(formulaNoSpaces) && unicode.IsDigit(rune(formulaNoSpaces[j])) {
					j++
				}
				i = j
			}
		} else if unicode.IsUpper(rune(formulaNoSpaces[i])) {
			i++
			if i < len(formulaNoSpaces) && unicode.IsLower(rune(formulaNoSpaces[i])) {
				i++
			}

			if i < len(formulaNoSpaces) && (formulaNoSpaces[i] == '-' || unicode.IsDigit(rune(formulaNoSpaces[i]))) {
				j := i
				if formulaNoSpaces[j] == '-' {
					j++
				}
				for j < len(formulaNoSpaces) && unicode.IsDigit(rune(formulaNoSpaces[j])) {
					j++
				}
				i = j
			}
		} else {
			return false
		}
	}

	return true
}

// validateGlycan statically validates a glycan string
func validateGlycan(glycan string) bool {
	glycanClean := strings.ReplaceAll(glycan, " ", "")

	monos := []string{
		"Hex", "HexNAc", "HexS", "HexP", "HexNAcS",
		"dHex", "NeuAc", "NeuGc", "Pen", "Fuc",
	}

	// Sort by length (longest first) to avoid partial matches
	sort.Slice(monos, func(i, j int) bool {
		return len(monos[i]) > len(monos[j])
	})

	// Build pattern for standard monosaccharides - count must start with 1-9
	monoPattern := "^("
	for i, mono := range monos {
		if i > 0 {
			monoPattern += "|"
		}
		monoPattern += regexp.QuoteMeta(mono)
	}
	monoPattern += `)((\([1-9]\d*\))|[1-9]\d*)?`

	// ProForma 2.1: Pattern for custom monosaccharides in curly braces
	// Format: {Formula} or {Formula:z+N} - count must start with 1-9
	customMonoPattern := `^\{([A-Za-z0-9]+)(:z[+-]\d+)?\}((\([1-9]\d*\))|[1-9]\d*)?`

	standardRe := regexp.MustCompile(monoPattern)
	customRe := regexp.MustCompile(customMonoPattern)

	i := 0
	for i < len(glycanClean) {
		// Try custom monosaccharide first
		customMatch := customRe.FindStringSubmatch(glycanClean[i:])
		if customMatch != nil {
			// Validate the formula part
			formula := customMatch[1]
			if !validateFormula(formula) {
				return false
			}

			monoLength := len(customMatch[0])
			hasCount := len(customMatch) > 3 && customMatch[3] != ""
			i += monoLength

			isAtEnd := i == len(glycanClean)
			if !hasCount && !isAtEnd {
				return false
			}
			continue
		}

		// Try standard monosaccharide
		standardMatch := standardRe.FindStringSubmatch(glycanClean[i:])
		if standardMatch != nil {
			monoLength := len(standardMatch[0])
			hasCount := len(standardMatch) > 2 && standardMatch[2] != ""
			i += monoLength

			isAtEnd := i == len(glycanClean)
			if !hasCount && !isAtEnd {
				return false
			}
			continue
		}

		// No match found
		return false
	}

	return i == len(glycanClean)
}
