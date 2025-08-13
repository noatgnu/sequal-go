package sequal

import (
	"encoding/json"
	"fmt"
	"sort"
)

// CalculateHashCode computes a hash code for a string (similar to TypeScript implementation)
func CalculateHashCode(input string) int32 {
	var hash int32 = 0
	for i := 0; i < len(input); i++ {
		char := int32(input[i])
		hash = (hash<<5 - hash) + char
		hash = hash & hash // Convert to 32-bit integer
	}
	return hash
}

// CountUniqueElements counts unique elements in a sequence of BaseBlock objects
func CountUniqueElements(seq []BaseBlock) map[string]int {
	elements := make(map[string]int)
	
	for _, item := range seq {
		// Count the base block value
		elements[item.GetValue()]++
		
		// Count modifications if the item is an amino acid with modifications
		if aminoAcid, ok := item.(*AminoAcid); ok {
			for _, mod := range aminoAcid.GetMods() {
				elements[mod.GetValue()]++
			}
		}
	}
	
	return elements
}

// VariablePositionPlacementGenerator generates all possible position combinations for modifications
func VariablePositionPlacementGenerator(positions []int) [][]int {
	if len(positions) == 0 {
		return [][]int{{}}
	}
	
	// Sort positions for consistent output
	sortedPositions := make([]int, len(positions))
	copy(sortedPositions, positions)
	sort.Ints(sortedPositions)
	
	var result [][]int
	n := len(sortedPositions)
	
	// Generate all possible combinations (2^n possibilities)
	for i := 0; i < (1 << n); i++ {
		var combination []int
		for j := 0; j < n; j++ {
			// Check if bit j is set in i
			if (i & (1 << j)) != 0 {
				combination = append(combination, sortedPositions[j])
			}
		}
		result = append(result, combination)
	}
	
	return result
}

// OrderedSerializePositionDict serializes a dictionary of positions with consistent ordering
func OrderedSerializePositionDict(positions map[int]interface{}) (string, error) {
	// Get all keys and sort them
	keys := make([]int, 0, len(positions))
	for k := range positions {
		keys = append(keys, k)
	}
	sort.Ints(keys)
	
	// Create ordered map
	sortedObj := make(map[string]interface{})
	for _, key := range keys {
		sortedObj[fmt.Sprintf("%d", key)] = positions[key]
	}
	
	jsonBytes, err := json.Marshal(sortedObj)
	if err != nil {
		return "", fmt.Errorf("could not serialize positions dictionary: %w", err)
	}
	
	return string(jsonBytes), nil
}

// SplitChimericProforma splits a chimeric ProForma string into individual peptidoforms
func SplitChimericProforma(proformaStr string) []string {
	var parts []string
	currentPartStart := 0
	bracketLevel := 0
	
	for i, char := range proformaStr {
		switch char {
		case '[', '{', '(':
			bracketLevel++
		case ']', '}', ')':
			if bracketLevel > 0 {
				bracketLevel--
			}
		case '+':
			if bracketLevel == 0 {
				// Found a separator '+' outside of any brackets
				part := proformaStr[currentPartStart:i]
				if len(part) > 0 {
					parts = append(parts, part)
				}
				currentPartStart = i + 1
			}
		}
	}
	
	// Add the last part of the string
	if currentPartStart < len(proformaStr) {
		part := proformaStr[currentPartStart:]
		if len(part) > 0 {
			parts = append(parts, part)
		}
	}
	
	return parts
}

// DeepCopyModifications creates deep copies of modification slices
func DeepCopyModifications(mods []*Modification) []*Modification {
	if mods == nil {
		return nil
	}
	
	result := make([]*Modification, len(mods))
	for i, mod := range mods {
		if mod != nil {
			// Create a copy of the modification
			result[i] = &Modification{}
			*result[i] = *mod
		}
	}
	
	return result
}

// IntPtr returns a pointer to an int value
func IntPtr(v int) *int {
	return &v
}

// Float64Ptr returns a pointer to a float64 value
func Float64Ptr(v float64) *float64 {
	return &v
}

// StringPtr returns a pointer to a string value
func StringPtr(v string) *string {
	return &v
}

// IntValue safely dereferences an int pointer, returning 0 if nil
func IntValue(ptr *int) int {
	if ptr == nil {
		return 0
	}
	return *ptr
}

// Float64Value safely dereferences a float64 pointer, returning 0.0 if nil
func Float64Value(ptr *float64) float64 {
	if ptr == nil {
		return 0.0
	}
	return *ptr
}

// StringValue safely dereferences a string pointer, returning empty string if nil
func StringValue(ptr *string) string {
	if ptr == nil {
		return ""
	}
	return *ptr
}