package sequal

import (
	"fmt"
	"regexp"
	"strconv"
	"strings"
)

// ProFormaParser handles parsing of ProForma 2.0 notation strings.
// It contains compiled regex patterns for efficient parsing of various ProForma elements.
type ProFormaParser struct {
	massShiftPattern    *regexp.Regexp
	crosslinkPattern    *regexp.Regexp
	crosslinkRefPattern *regexp.Regexp
	branchPattern       *regexp.Regexp
	branchRefPattern    *regexp.Regexp
}

// NewProFormaParser creates a new ProFormaParser with pre-compiled regex patterns
// for parsing mass shifts, crosslinks, and branches.
func NewProFormaParser() *ProFormaParser {
	return &ProFormaParser{
		massShiftPattern:    regexp.MustCompile(`^[+-]\d+(\.\d+)?$`),
		crosslinkPattern:    regexp.MustCompile(`^([^#]+)#(XL[A-Za-z0-9]+)$`),
		crosslinkRefPattern: regexp.MustCompile(`^#(XL[A-Za-z0-9]+)$`),
		branchPattern:       regexp.MustCompile(`^([^#]+)#BRANCH$`),
		branchRefPattern:    regexp.MustCompile(`^#BRANCH$`),
	}
}

// ParseProFormaResult contains all parsed components from a ProForma string
// including the base sequence, modifications, global modifications, sequence ambiguities,
// charge state, and ionic species information.
type ParseProFormaResult struct {
	BaseSequence        string
	Modifications       map[string][]*Modification
	GlobalMods          []*GlobalModification
	SequenceAmbiguities []*SequenceAmbiguity
	Charge              *int
	IonicSpecies        *string
}

// ParseProForma parses a ProForma string and returns its basic components.
// This is a convenience function that creates a parser and calls Parse.
func ParseProForma(proformaStr string) (string, map[string][]*Modification, []*GlobalModification, []*SequenceAmbiguity, []*int, error) {
	parser := NewProFormaParser()
	return parser.Parse(proformaStr)
}

// ParseProFormaDetailed parses a ProForma string and returns a structured result
// containing all parsed components including charge and ionic species information.
func ParseProFormaDetailed(proformaStr string) (*ParseProFormaResult, error) {
	parser := NewProFormaParser()
	baseSeq, mods, globalMods, seqAmbig, chargeInfo, err := parser.Parse(proformaStr)
	if err != nil {
		return nil, err
	}
	
	result := &ParseProFormaResult{
		BaseSequence:        baseSeq,
		Modifications:       mods,
		GlobalMods:          globalMods,
		SequenceAmbiguities: seqAmbig,
	}
	
	if len(chargeInfo) > 0 {
		result.Charge = chargeInfo[0]
	}
	
	if strings.Contains(proformaStr, "/") {
		chargeInfoResult, err := parser.parseChargeInfo(proformaStr)
		if err == nil && len(chargeInfoResult) > 2 {
			if species, ok := chargeInfoResult[2].(*string); ok {
				result.IonicSpecies = species
			}
		}
	}
	
	return result, nil
}

// Parse parses a ProForma string into its constituent parts and returns the base sequence,
// modifications map, global modifications, sequence ambiguities, and charge information.
// This is the main parsing method that handles all ProForma 2.0 notation elements.
func (p *ProFormaParser) Parse(proformaStr string) (string, map[string][]*Modification, []*GlobalModification, []*SequenceAmbiguity, []*int, error) {
	baseSequence := ""
	modifications := make(map[string][]*Modification)
	globalMods := make([]*GlobalModification, 0)
	sequenceAmbiguities := make([]*SequenceAmbiguity, 0)

	getModsAtPosition := func(pos int) []*Modification {
		posStr := strconv.Itoa(pos)
		if modifications[posStr] == nil {
			modifications[posStr] = make([]*Modification, 0)
		}
		return modifications[posStr]
	}

	setModsAtPosition := func(pos int, mods []*Modification) {
		posStr := strconv.Itoa(pos)
		modifications[posStr] = mods
	}

	// Parse global modifications
	for strings.HasPrefix(proformaStr, "<") {
		endBracket := strings.Index(proformaStr, ">")
		if endBracket == -1 {
			return "", nil, nil, nil, nil, fmt.Errorf("unclosed global modification angle bracket")
		}

		globalModStr := proformaStr[1:endBracket]
		proformaStr = proformaStr[endBracket+1:]

		if strings.Contains(globalModStr, "@") {
			// Fixed protein modification
			parts := strings.Split(globalModStr, "@")
			if len(parts) != 2 {
				return "", nil, nil, nil, nil, fmt.Errorf("invalid global modification format")
			}

			modPart, targets := parts[0], parts[1]
			modValue := modPart

			if strings.HasPrefix(modPart, "[") && strings.HasSuffix(modPart, "]") {
				modValue = modPart[1 : len(modPart)-1]
			}

			targetResidues := strings.Split(targets, ",")
			globalMods = append(globalMods, NewGlobalModification(modValue, targetResidues, "fixed"))
		} else {
			// Isotope labeling
			globalMods = append(globalMods, NewGlobalModification(globalModStr, nil, "isotope"))
		}
	}

	// Handle unknown position modifications
	if strings.Contains(proformaStr, "?") {
		i := 0
		var unknownPosMods []string
		proformaRunes := []rune(proformaStr)

		for i < len(proformaRunes) {
			if proformaRunes[i] != '[' {
				if len(unknownPosMods) > 0 && i < len(proformaRunes) && proformaRunes[i] == '?' {
					for _, modStr := range unknownPosMods {
						mod := p.createModification(modStr, map[string]interface{}{"isUnknownPosition": true})
						currentMods := getModsAtPosition(-4)
						currentMods = append(currentMods, mod)
						setModsAtPosition(-4, currentMods)
					}
					i++
				}
				unknownPosMods = nil
				break
			}

			bracketCount := 1
			j := i + 1
			for j < len(proformaRunes) && bracketCount > 0 {
				if proformaRunes[j] == '[' {
					bracketCount++
				} else if proformaRunes[j] == ']' {
					bracketCount--
				}
				j++
			}

			if bracketCount > 0 {
				return "", nil, nil, nil, nil, fmt.Errorf("unclosed bracket at position %d", i)
			}

			modStr := string(proformaRunes[i+1 : j-1])

			count := 1
			if j < len(proformaRunes) && proformaRunes[j] == '^' {
				j++
				numStart := j
				for j < len(proformaRunes) && proformaRunes[j] >= '0' && proformaRunes[j] <= '9' {
					j++
				}
				if j > numStart {
					var err error
					count, err = strconv.Atoi(string(proformaRunes[numStart:j]))
					if err != nil {
						count = 1
					}
				}
			}

			for k := 0; k < count; k++ {
				unknownPosMods = append(unknownPosMods, modStr)
			}
			i = j
		}
		proformaStr = string(proformaRunes[i:])
	}

	// Parse labile modifications
	i := 0
	for i < len(proformaStr) && proformaStr[i] == '{' {
		j := strings.Index(proformaStr[i:], "}")
		if j == -1 {
			return "", nil, nil, nil, nil, fmt.Errorf("unclosed curly brace at position %d", i)
		}
		j += i

		modStr := proformaStr[i+1 : j]
		if !strings.HasPrefix(modStr, "Glycan:") {
			return "", nil, nil, nil, nil, fmt.Errorf("labile modification must start with 'Glycan:', found: %s", modStr)
		}

		mod := p.createModification(modStr, map[string]interface{}{"isLabile": true})
		currentMods := getModsAtPosition(-3)
		currentMods = append(currentMods, mod)
		setModsAtPosition(-3, currentMods)
		i = j + 1
	}

	proformaStr = proformaStr[i:]

	// Parse N-terminal modifications
	if strings.HasPrefix(proformaStr, "[") {
		bracketLevel := 0
		terminatorPos := -1

		for i, char := range proformaStr {
			switch char {
			case '[':
				bracketLevel++
			case ']':
				bracketLevel--
			case '-':
				if bracketLevel == 0 {
					terminatorPos = i
					break
				}
			}
			if terminatorPos != -1 {
				break
			}
		}

		if terminatorPos != -1 {
			nTerminalPart := proformaStr[:terminatorPos]
			proformaStr = proformaStr[terminatorPos+1:]

			currentPos := 0
			for currentPos < len(nTerminalPart) {
				if nTerminalPart[currentPos] == '[' {
					bracketDepth := 1
					endPos := currentPos + 1

					for endPos < len(nTerminalPart) && bracketDepth > 0 {
						if nTerminalPart[endPos] == '[' {
							bracketDepth++
						}
						if nTerminalPart[endPos] == ']' {
							bracketDepth--
						}
						endPos++
					}

					if bracketDepth == 0 {
						modString := nTerminalPart[currentPos+1 : endPos-1]
						nTermMod := p.createModification(modString, map[string]interface{}{"isTerminal": true})
						currentMods := getModsAtPosition(-1)
						currentMods = append(currentMods, nTermMod)
						setModsAtPosition(-1, currentMods)
					}

					currentPos = endPos
				} else {
					currentPos++
				}
			}
		}
	}

	// Parse charge information
	chargeInfo, err := p.parseChargeInfo(proformaStr)
	if err != nil {
		return "", nil, nil, nil, nil, err
	}
	proformaStr = chargeInfo[0].(string)

	// Parse C-terminal modifications
	if strings.Contains(proformaStr, "-") {
		bracketLevel := 0
		terminatorPos := -1

		// Scan from right to left
		proformaRunes := []rune(proformaStr)
		for i := len(proformaRunes) - 1; i >= 0; i-- {
			switch proformaRunes[i] {
			case ']':
				bracketLevel++
			case '[':
				bracketLevel--
			case '-':
				if bracketLevel == 0 {
					terminatorPos = i
					break
				}
			}
			if terminatorPos != -1 {
				break
			}
		}

		if terminatorPos != -1 {
			cTerminalPart := string(proformaRunes[terminatorPos+1:])
			proformaStr = string(proformaRunes[:terminatorPos])

			currentPos := 0
			for currentPos < len(cTerminalPart) {
				if cTerminalPart[currentPos] == '[' {
					bracketDepth := 1
					endPos := currentPos + 1

					for endPos < len(cTerminalPart) && bracketDepth > 0 {
						if cTerminalPart[endPos] == '[' {
							bracketDepth++
						}
						if cTerminalPart[endPos] == ']' {
							bracketDepth--
						}
						endPos++
					}

					if bracketDepth == 0 {
						modString := cTerminalPart[currentPos+1 : endPos-1]
						cTermMod := p.createModification(modString, map[string]interface{}{"isTerminal": true})
						currentMods := getModsAtPosition(-2)
						currentMods = append(currentMods, cTermMod)
						setModsAtPosition(-2, currentMods)
					}

					currentPos = endPos
				} else {
					currentPos++
				}
			}
		}
	}

	// Parse main sequence
	i = 0
	nextModIsGap := false
	var rangeStack []int

	for i < len(proformaStr) {
		char := proformaStr[i]

		if i+1 < len(proformaStr) && proformaStr[i:i+2] == "(?" {
			closingParen := strings.Index(proformaStr[i+2:], ")")
			if closingParen == -1 {
				return "", nil, nil, nil, nil, fmt.Errorf("unclosed sequence ambiguity parenthesis")
			}
			closingParen += i + 2

			ambiguousSeq := proformaStr[i+2 : closingParen]
			sequenceAmbiguities = append(sequenceAmbiguities, NewSequenceAmbiguity(ambiguousSeq, 0))

			i = closingParen + 1
			continue
		}

		switch char {
		case '(':
			rangeStack = append(rangeStack, len(baseSequence))
			i++
			continue

		case ')':
			if len(rangeStack) == 0 {
				return "", nil, nil, nil, nil, fmt.Errorf("unmatched closing parenthesis")
			}

			rangeStart := rangeStack[len(rangeStack)-1]
			rangeStack = rangeStack[:len(rangeStack)-1]
			rangeEnd := len(baseSequence) - 1

			// Look for modification after the range
			j := i + 1
			for j < len(proformaStr) && proformaStr[j] == '[' {
				modStart := j
				bracketCount := 1
				j++

				for j < len(proformaStr) && bracketCount > 0 {
					if proformaStr[j] == '[' {
						bracketCount++
					} else if proformaStr[j] == ']' {
						bracketCount--
					}
					j++
				}

				if bracketCount == 0 {
					modStr := proformaStr[modStart+1 : j-1]
					mod := p.createModification(modStr, map[string]interface{}{
						"inRange":    true,
						"rangeStart": rangeStart,
						"rangeEnd":   rangeEnd,
					})

					for pos := rangeStart; pos <= rangeEnd; pos++ {
						currentMods := getModsAtPosition(pos)
						currentMods = append(currentMods, mod)
						setModsAtPosition(pos, currentMods)
					}
				}
			}
			i = j

		case '[':
			bracketCount := 1
			j := i + 1
			for j < len(proformaStr) && bracketCount > 0 {
				if proformaStr[j] == '[' {
					bracketCount++
				} else if proformaStr[j] == ']' {
					bracketCount--
				}
				j++
			}

			if bracketCount > 0 {
				return "", nil, nil, nil, nil, fmt.Errorf("unclosed square bracket at position %d", i)
			}

			modStr := proformaStr[i+1 : j-1]
			var mod *Modification

			if nextModIsGap {
				mod = p.createModification(modStr, map[string]interface{}{"isGap": true})
				nextModIsGap = false
			} else if p.crosslinkRefPattern.MatchString(modStr) {
				mod = p.createModification(modStr, map[string]interface{}{"isCrosslinkRef": true})
			} else if p.branchRefPattern.MatchString(modStr) {
				mod = p.createModification(modStr, map[string]interface{}{"isBranchRef": true})
			} else if matches := p.crosslinkPattern.FindStringSubmatch(modStr); matches != nil {
				mod = p.createModification(modStr, map[string]interface{}{"crosslinkId": matches[2]})
			} else if matches := p.branchPattern.FindStringSubmatch(modStr); matches != nil {
				mod = p.createModification(modStr, map[string]interface{}{"isBranch": true})
			} else {
				mod = p.createModification(modStr, nil)
			}

			if len(baseSequence) > 0 {
				currentMods := getModsAtPosition(len(baseSequence) - 1)
				currentMods = append(currentMods, mod)
				setModsAtPosition(len(baseSequence)-1, currentMods)
			}

			i = j

		case '{':
			j := strings.Index(proformaStr[i:], "}")
			if j == -1 {
				return "", nil, nil, nil, nil, fmt.Errorf("unclosed curly brace at position %d", i)
			}
			j += i

			modStr := proformaStr[i+1 : j]
			mod := p.createModification(modStr, map[string]interface{}{"isAmbiguous": true})

			if len(baseSequence) > 0 {
				currentMods := getModsAtPosition(len(baseSequence) - 1)
				currentMods = append(currentMods, mod)
				setModsAtPosition(len(baseSequence)-1, currentMods)
			}

			i = j + 1

		default:
			baseSequence += string(char)
			isGap := char == 'X' && i+1 < len(proformaStr) && proformaStr[i+1] == '['
			if isGap {
				nextModIsGap = true
			}
			i++
		}
	}

	if len(rangeStack) > 0 {
		return "", nil, nil, nil, nil, fmt.Errorf("unclosed parenthesis")
	}

	var chargeInfoResult []*int
	if len(chargeInfo) > 1 {
		if charge, ok := chargeInfo[1].(*int); ok && charge != nil {
			chargeInfoResult = append(chargeInfoResult, charge)
		} else {
			chargeInfoResult = append(chargeInfoResult, nil)
		}
	}

	return baseSequence, modifications, globalMods, sequenceAmbiguities, chargeInfoResult, nil
}

// createModification creates a Modification instance with the specified options.
// The options map contains various boolean flags and values that control the modification type.
func (p *ProFormaParser) createModification(modStr string, options map[string]interface{}) *Modification {
	isTerminal := false
	isAmbiguous := false
	isLabile := false
	isUnknownPosition := false
	var crosslinkId *string
	isCrosslinkRef := false
	isBranch := false
	isBranchRef := false
	isGap := false
	inRange := false
	var rangeStart, rangeEnd *int

	// Extract options
	if options != nil {
		if v, ok := options["isTerminal"].(bool); ok {
			isTerminal = v
		}
		if v, ok := options["isAmbiguous"].(bool); ok {
			isAmbiguous = v
		}
		if v, ok := options["isLabile"].(bool); ok {
			isLabile = v
		}
		if v, ok := options["isUnknownPosition"].(bool); ok {
			isUnknownPosition = v
		}
		if v, ok := options["crosslinkId"].(string); ok {
			crosslinkId = &v
		}
		if v, ok := options["isCrosslinkRef"].(bool); ok {
			isCrosslinkRef = v
		}
		if v, ok := options["isBranch"].(bool); ok {
			isBranch = v
		}
		if v, ok := options["isBranchRef"].(bool); ok {
			isBranchRef = v
		}
		if v, ok := options["isGap"].(bool); ok {
			isGap = v
		}
		if v, ok := options["inRange"].(bool); ok {
			inRange = v
		}
		if v, ok := options["rangeStart"].(int); ok {
			rangeStart = &v
		}
		if v, ok := options["rangeEnd"].(int); ok {
			rangeEnd = &v
		}
	}

	modValue := NewModificationValue(modStr, nil)
	modType := "static"

	if isTerminal {
		modType = "terminal"
	} else if isAmbiguous {
		modType = "ambiguous"
	} else if isLabile {
		modType = "labile"
	} else if isUnknownPosition {
		modType = "unknown_position"
	} else if crosslinkId != nil || isCrosslinkRef {
		modType = "crosslink"
	} else if isBranch || isBranchRef {
		modType = "branch"
	} else if isGap {
		modType = "gap"
	}

	// Handle mass shifts
	if p.massShiftPattern.MatchString(modStr) && !strings.Contains(modStr, "#") {
		massValue, _ := strconv.ParseFloat(modStr, 64)
		modValueForMassShift := NewModificationValue("Mass:"+modStr, &massValue)
		if isGap {
			return NewModification(modStr, nil, nil, nil, "gap", false, 0, massValue, false,
				nil, false, false, false, nil, false, inRange, rangeStart, rangeEnd, nil, modValueForMassShift)
		} else if inRange {
			return NewModification(modStr, nil, nil, nil, "variable", false, 0, massValue, false,
				nil, false, false, false, nil, false, true, rangeStart, rangeEnd, nil, modValueForMassShift)
		}
		return NewModification("Mass:"+modStr, nil, nil, nil, "static", false, 0, massValue, false,
			nil, false, false, false, nil, false, inRange, rangeStart, rangeEnd, nil, modValueForMassShift)
	}

	// Handle ambiguity patterns
	ambiguityPattern := regexp.MustCompile(`(.+?)#([A-Za-z0-9]+)(?:\(([0-9.]+)\))?$`)
	ambiguityRefPattern := regexp.MustCompile(`#([A-Za-z0-9]+)(?:\(([0-9.]+)\))?$`)

	if strings.Contains(modStr, "#") && !isCrosslinkRef && !isBranch && !isBranchRef && crosslinkId == nil {
		if matches := ambiguityPattern.FindStringSubmatch(modStr); matches != nil && !strings.HasPrefix(matches[2], "XL") {
			modStr = matches[1]
			ambiguityGroup := matches[2]
			var localizationScore *float64
			if len(matches) > 3 && matches[3] != "" {
				if score, err := strconv.ParseFloat(matches[3], 64); err == nil {
					localizationScore = &score
				}
			}
			return NewModification(modStr, nil, nil, nil, "ambiguous", false, 0, 0.0, false,
				nil, false, false, false, &ambiguityGroup, false, inRange, rangeStart, rangeEnd, localizationScore, modValue)
		} else if matches := ambiguityRefPattern.FindStringSubmatch(modStr); matches != nil && !strings.HasPrefix(matches[1], "XL") {
			ambiguityGroup := matches[1]
			var localizationScore *float64
			if len(matches) > 2 && matches[2] != "" {
				if score, err := strconv.ParseFloat(matches[2], 64); err == nil {
					localizationScore = &score
				}
			}
			return NewModification("", nil, nil, nil, "ambiguous", false, 0, 0.0, false,
				nil, false, false, false, &ambiguityGroup, true, inRange, rangeStart, rangeEnd, localizationScore, modValue)
		}
	}

	// Create the modification with appropriate attributes
	return NewModification(modStr, nil, nil, nil, modType, isLabile, 0, 0.0, false,
		crosslinkId, isCrosslinkRef, isBranchRef, isBranch, nil, false, inRange, rangeStart, rangeEnd, nil, modValue)
}

// parseChargeInfo parses charge information from a ProForma string.
// Returns the modified string (without charge info), charge value, and ionic species.
func (p *ProFormaParser) parseChargeInfo(proformaStr string) ([]interface{}, error) {
	if !strings.Contains(proformaStr, "/") {
		return []interface{}{proformaStr, nil, nil}, nil
	}

	chargePos := -1
	bracketLevel := 0
	for i, char := range proformaStr {
		switch char {
		case '[', '(':
			bracketLevel++
		case ']', ')':
			bracketLevel--
		case '/':
			if bracketLevel == 0 {
				chargePos = i
				break
			}
		}
		if chargePos != -1 {
			break
		}
	}

	if chargePos == -1 {
		return []interface{}{proformaStr, nil, nil}, nil
	}

	beforeCharge := proformaStr[:chargePos]
	afterCharge := proformaStr[chargePos+1:]

	i := 0
	sign := 1

	// Handle negative sign
	if i < len(afterCharge) && afterCharge[i] == '-' {
		sign = -1
		i++
	}

	// Parse digits
	startDigit := i
	for i < len(afterCharge) && afterCharge[i] >= '0' && afterCharge[i] <= '9' {
		i++
	}

	if startDigit == i { // No digits found
		return []interface{}{proformaStr, nil, nil}, nil
	}

	chargeValue, err := strconv.Atoi(afterCharge[startDigit:i])
	if err != nil {
		return []interface{}{proformaStr, nil, nil}, nil
	}
	chargeValue *= sign

	// Check for ionic species in square brackets
	remaining := afterCharge[i:]
	var ionicSpecies *string

	if len(remaining) > 0 && remaining[0] == '[' {
		// Find the matching closing bracket
		bracketLevel := 1
		endPos := 0

		for j := 1; j < len(remaining); j++ {
			if remaining[j] == '[' {
				bracketLevel++
			} else if remaining[j] == ']' {
				bracketLevel--
			}

			if bracketLevel == 0 {
				endPos = j
				break
			}
		}

		if endPos > 0 {
			species := remaining[1:endPos]
			ionicSpecies = &species
			remaining = remaining[endPos+1:]
		}
	}

	// Reconstruct the string without charge information
	resultStr := beforeCharge
	if len(remaining) > 0 {
		resultStr += remaining
	}

	return []interface{}{resultStr, &chargeValue, ionicSpecies}, nil
}