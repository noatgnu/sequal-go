package sequal

import (
	"fmt"
	"regexp"
	"strconv"
	"strings"
)

// Sequence represents a peptide sequence with modifications and supports ProForma notation.
// It can handle single sequences, multi-chain complexes, and chimeric peptidoforms.
type Sequence struct {
	seq                 []*AminoAcid
	chains              []*Sequence
	isMultiChain        bool
	mods                map[int][]*Modification
	globalMods          []*GlobalModification
	sequenceAmbiguities []*SequenceAmbiguity
	seqLength           int
	charge              *int
	ionicSpecies        *string
	isChimeric          bool
	peptidoforms        []*Sequence
}

// NewSequence creates a new Sequence instance with the specified parameters.
// The seq parameter can be a string or other sequence representation.
// If parse is true, the sequence will be parsed according to modPosition ("left" or "right").
func NewSequence(seq interface{}, mods map[int][]*Modification, parse bool, 
	modPosition string, chains []*Sequence, globalMods []*GlobalModification, 
	sequenceAmbiguities []*SequenceAmbiguity, charge *int, ionicSpecies *string) *Sequence {
	
	s := &Sequence{
		seq:                 make([]*AminoAcid, 0),
		chains:              chains,
		isMultiChain:        false,
		mods:                make(map[int][]*Modification),
		globalMods:          globalMods,
		sequenceAmbiguities: sequenceAmbiguities,
		charge:              charge,
		ionicSpecies:        ionicSpecies,
		isChimeric:          false,
		peptidoforms:        make([]*Sequence, 0),
	}

	if mods != nil {
		for pos, modList := range mods {
			s.mods[pos] = DeepCopyModifications(modList)
		}
	}

	if parse {
		s.parseSequence(seq, modPosition)
	}

	s.seqLength = len(s.seq)
	return s
}

// FromProforma creates a Sequence object from a ProForma notation string.
// Supports all ProForma 2.0 features including multi-chain sequences (//), 
// chimeric sequences (+), and all modification types.
func FromProforma(proformaStr string) (*Sequence, error) {
	if strings.Contains(proformaStr, "//") {
		chains := strings.Split(proformaStr, "//")
		mainSeq, err := FromProforma(chains[0])
		if err != nil {
			return nil, err
		}
		mainSeq.isMultiChain = true
		mainSeq.chains = []*Sequence{mainSeq}

		for i := 1; i < len(chains); i++ {
			chain, err := FromProforma(chains[i])
			if err != nil {
				return nil, err
			}
			mainSeq.chains = append(mainSeq.chains, chain)
		}

		return mainSeq, nil
	}

	peptidoforms := SplitChimericProforma(proformaStr)
	if len(peptidoforms) > 1 {
		mainSeq, err := FromProforma(peptidoforms[0])
		if err != nil {
			return nil, err
		}
		mainSeq.isChimeric = true
		mainSeq.peptidoforms = []*Sequence{mainSeq}

		for i := 1; i < len(peptidoforms); i++ {
			peptidoform, err := FromProforma(peptidoforms[i])
			if err != nil {
				return nil, err
			}
			peptidoform.isChimeric = true
			mainSeq.peptidoforms = append(mainSeq.peptidoforms, peptidoform)
		}

		return mainSeq, nil
	}
	result, err := ParseProFormaDetailed(proformaStr)
	if err != nil {
		return nil, err
	}

	baseSequence := result.BaseSequence
	modifications := result.Modifications
	globalMods := result.GlobalMods
	sequenceAmbiguities := result.SequenceAmbiguities
	charge := result.Charge
	species := result.IonicSpecies

	seq := NewSequence(
		baseSequence,
		make(map[int][]*Modification),
		true,
		"right",
		[]*Sequence{},
		globalMods,
		sequenceAmbiguities,
		charge,
		species,
	)

	if charge != nil {
		seq.isChimeric = true
	}

	for posStr, mods := range modifications {
		pos, _ := strconv.Atoi(posStr)
		for _, mod := range mods {
			switch pos {
			case -1, -2, -3, -4:
				if seq.mods[pos] == nil {
					seq.mods[pos] = make([]*Modification, 0)
				}
				seq.mods[pos] = append(seq.mods[pos], mod)
			default:
				if pos >= 0 && pos < len(seq.seq) {
					seq.seq[pos].AddModification(mod)
				}
			}
		}
	}

	seq.peptidoforms = []*Sequence{seq}
	return seq, nil
}

// parseSequence parses the input sequence into a list of AminoAcid objects.
// The modPosition parameter determines whether modifications are applied to the
// left ("left") or right ("right") of the amino acid they modify.
func (s *Sequence) parseSequence(seq interface{}, modPosition string) error {
	if modPosition != "left" && modPosition != "right" {
		return fmt.Errorf("modPosition must be either 'left' or 'right'")
	}

	var currentMod []*Modification
	currentPosition := 0

	seqStr := ""
	switch v := seq.(type) {
	case string:
		seqStr = v
	case []interface{}:
		for _, item := range v {
			if str, ok := item.(string); ok {
				seqStr += str
			}
		}
	default:
		return fmt.Errorf("unsupported sequence type")
	}

	for _, block := range s.sequenceIterator(seqStr) {
		if !block.IsMod {
			if modPosition == "left" {
				aa, err := NewAminoAcid(block.Value, &currentPosition, nil)
				if err != nil {
					return err
				}
				
				for _, mod := range currentMod {
					aa.AddModification(mod)
				}
				s.seq = append(s.seq, aa)
				currentMod = nil
			} else {
				if len(currentMod) > 0 && currentPosition > 0 {
					for _, mod := range currentMod {
						s.seq[currentPosition-1].AddModification(mod)
					}
				}

				aa, err := NewAminoAcid(block.Value, &currentPosition, nil)
				if err != nil {
					return err
				}

				if mods, exists := s.mods[currentPosition]; exists {
					for _, mod := range mods {
						aa.AddModification(mod)
					}
				}

				s.seq = append(s.seq, aa)
				currentMod = nil
			}
			currentPosition++
		} else {
			// Handle modification
			if len(s.mods) == 0 {
				modValue := s.extractModValue(block.Value)
				mod := NewModification(modValue, nil, nil, nil, "static", false, 0, 0.0, false,
					nil, false, false, false, nil, false, false, nil, nil, nil, nil)

				if modPosition == "right" && currentPosition > 0 {
					s.seq[currentPosition-1].AddModification(mod)
				} else {
					currentMod = append(currentMod, mod)
				}
			}
		}
	}

	return nil
}

// SequenceBlock represents a parsed block from sequence iteration
type SequenceBlock struct {
	Value string
	IsMod bool
}

// sequenceIterator iterates through sequence elements, identifying blocks and modifications
func (s *Sequence) sequenceIterator(seq string) []SequenceBlock {
	var result []SequenceBlock
	modOpen := 0
	block := ""
	isMod := false

	modEnclosureStart := map[rune]bool{'(': true, '[': true, '{': true}
	modEnclosureEnd := map[rune]bool{')': true, ']': true, '}': true}

	for _, char := range seq {
		if modEnclosureStart[char] {
			isMod = true
			modOpen++
		} else if modEnclosureEnd[char] {
			modOpen--
		}
		block += string(char)

		if modOpen == 0 && block != "" {
			result = append(result, SequenceBlock{Value: block, IsMod: isMod})
			isMod = false
			block = ""
		}
	}

	return result
}

// extractModValue extracts modification value from a string
func (s *Sequence) extractModValue(modStr string) string {
	if len(modStr) >= 2 {
		start := modStr[0]
		end := modStr[len(modStr)-1]
		if (start == '(' && end == ')') || (start == '[' && end == ']') || (start == '{' && end == '}') {
			return modStr[1 : len(modStr)-1]
		}
	}
	return modStr
}

// ToProforma converts the sequence to ProForma format
func (s *Sequence) ToProforma() string {
	if s.isMultiChain {
		chains := make([]string, len(s.chains))
		for i, chain := range s.chains {
			chains[i] = s.chainToProforma(chain)
		}
		return strings.Join(chains, "//")
	} else if s.isChimeric && len(s.peptidoforms) > 0 {
		peptidoforms := make([]string, len(s.peptidoforms))
		for i, pep := range s.peptidoforms {
			peptidoforms[i] = s.chainToProforma(pep)
		}
		return strings.Join(peptidoforms, "+")
	}
	return s.chainToProforma(s)
}

// chainToProforma converts a chain to ProForma format
func (s *Sequence) chainToProforma(chain *Sequence) string {
	result := ""

	// Add global modifications
	for _, mod := range s.globalMods {
		result += mod.ToProforma()
	}

	// Handle unknown position modifications (-4)
	if unknownMods, exists := chain.mods[-4]; exists {
		unknownModsByValue := make(map[string]int)
		for _, mod := range unknownMods {
			modProforma := mod.ToProforma()
			unknownModsByValue[modProforma]++
		}

		for modValue, count := range unknownModsByValue {
			if count > 1 {
				result += fmt.Sprintf("[%s]^%d?", modValue, count)
			} else {
				result += fmt.Sprintf("[%s]?", modValue)
			}
		}
	}

	// Handle labile modifications (-3)
	if labileMods, exists := chain.mods[-3]; exists {
		for _, mod := range labileMods {
			if mod.GetModType() == "labile" {
				result += fmt.Sprintf("{%s}", mod.ToProforma())
			}
		}
	}

	// Handle N-terminal modifications (-1)
	if nTermMods, exists := chain.mods[-1]; exists {
		nModStr := ""
		for _, mod := range nTermMods {
			nModStr += fmt.Sprintf("[%s]", mod.ToProforma())
		}
		if nModStr != "" {
			result += nModStr + "-"
		}
	}

	// Process each amino acid in the sequence
	for _, aa := range chain.seq {
		// Add amino acid value
		result += aa.GetValue()

		// Add modifications for this position
		mods := aa.GetMods()
		if len(mods) > 0 {
			for _, mod := range mods {
				modStr := mod.ToProforma()
				if mod.GetModType() == "ambiguous" && !mod.HasAmbiguity() {
					// Use curly braces for ambiguous modifications without ambiguity groups
					result += fmt.Sprintf("{%s}", modStr)
				} else {
					// Use square brackets for all other modifications
					result += fmt.Sprintf("[%s]", modStr)
				}
			}
		}
	}

	// Handle C-terminal modifications (-2)
	if cTermMods, exists := chain.mods[-2]; exists {
		cModStr := ""
		for _, mod := range cTermMods {
			cModStr += fmt.Sprintf("[%s]", mod.ToProforma())
		}
		if cModStr != "" {
			result += "-" + cModStr
		}
	}

	// Add charge information
	if chain.charge != nil {
		result += fmt.Sprintf("/%d", *chain.charge)
		if chain.ionicSpecies != nil {
			result += fmt.Sprintf("[%s]", *chain.ionicSpecies)
		}
	}

	return result
}

// ToStrippedString returns the sequence as a string without any modification annotations
func (s *Sequence) ToStrippedString() string {
	result := ""
	for _, aa := range s.seq {
		result += aa.GetValue()
	}
	return result
}

// GetLength returns the length of the sequence
func (s *Sequence) GetLength() int {
	return s.seqLength
}

// String returns a string representation of the sequence
func (s *Sequence) String() string {
	result := ""
	for _, aa := range s.seq {
		result += aa.String()
	}
	return result
}

// GetItem returns an amino acid at the specified position or a slice
func (s *Sequence) GetItem(key interface{}) interface{} {
	switch k := key.(type) {
	case int:
		if k >= 0 && k < len(s.seq) {
			return s.seq[k]
		}
		return nil
	case []int:
		if len(k) == 2 {
			start, end := k[0], k[1]
			if start >= 0 && end <= len(s.seq) && start <= end {
				newSeq := NewSequence("", nil, false, "right", nil, nil, nil, nil, nil)
				newSeq.seq = s.seq[start:end]
				newSeq.seqLength = len(newSeq.seq)
				return newSeq
			}
		}
		return nil
	}
	return nil
}

// Equal checks if two sequences are equal
func (s *Sequence) Equal(other *Sequence) bool {
	if other == nil {
		return false
	}
	if s.seqLength != other.seqLength {
		return false
	}
	for i := 0; i < s.seqLength; i++ {
		if !s.seq[i].Equal(other.seq[i]) {
			return false
		}
	}
	return true
}

// AddModifications adds modifications to residues at specified positions
func (s *Sequence) AddModifications(modDict map[int][]*Modification) {
	for _, aa := range s.seq {
		if aa.GetPosition() != nil {
			pos := *aa.GetPosition()
			if mods, exists := modDict[pos]; exists {
				for _, mod := range mods {
					aa.AddModification(mod)
				}
			}
		}
	}
}

// FindWithRegex finds positions in the sequence that match a given regex motif
func (s *Sequence) FindWithRegex(motif string, ignore []bool) ([][]int, error) {
	pattern, err := regexp.Compile(motif)
	if err != nil {
		return nil, err
	}

	var seqStr string
	if ignore != nil && len(ignore) == len(s.seq) {
		// Build string excluding ignored positions
		for i, aa := range s.seq {
			if !ignore[i] {
				seqStr += aa.GetValue()
			}
		}
	} else {
		seqStr = s.ToStrippedString()
	}

	var results [][]int
	matches := pattern.FindAllStringSubmatchIndex(seqStr, -1)

	for _, match := range matches {
		if len(match) > 2 {
			// Has capture groups
			for i := 1; i < len(match)/2; i++ {
				start := match[i*2]
				end := match[i*2+1]
				if start >= 0 && end >= 0 {
					results = append(results, []int{start, end})
				}
			}
		} else {
			results = append(results, []int{match[0], match[1]})
		}
	}

	return results, nil
}

// Gaps identifies gaps in the sequence
func (s *Sequence) Gaps() []bool {
	gaps := make([]bool, len(s.seq))
	for i, aa := range s.seq {
		gaps[i] = aa.GetValue() == "-"
	}
	return gaps
}

// Count counts occurrences of a character in a range
func (s *Sequence) Count(char string, start, end int) int {
	subStr := s.ToStrippedString()
	if end > len(subStr) {
		end = len(subStr)
	}
	if start < 0 {
		start = 0
	}
	if start >= end {
		return 0
	}

	subStr = subStr[start:end]
	return strings.Count(subStr, char)
}

// ToMap converts the sequence to a map representation
func (s *Sequence) ToMap() map[string]interface{} {
	modsByPosition := make(map[string]interface{})

	for i, aa := range s.seq {
		mods := aa.GetMods()
		if len(mods) > 0 {
			modMaps := make([]map[string]interface{}, len(mods))
			for j, mod := range mods {
				modMaps[j] = mod.ToMap()
			}
			modsByPosition[strconv.Itoa(i)] = modMaps
		}
	}

	return map[string]interface{}{
		"sequence":      s.ToStrippedString(),
		"modifications": modsByPosition,
	}
}

// GetSeq returns the amino acid sequence
func (s *Sequence) GetSeq() []*AminoAcid {
	return s.seq
}

// GetMods returns the modifications map
func (s *Sequence) GetMods() map[int][]*Modification {
	return s.mods
}

// GetGlobalMods returns the global modifications
func (s *Sequence) GetGlobalMods() []*GlobalModification {
	return s.globalMods
}

// GetSequenceAmbiguities returns the sequence ambiguities
func (s *Sequence) GetSequenceAmbiguities() []*SequenceAmbiguity {
	return s.sequenceAmbiguities
}

// GetCharge returns the charge state
func (s *Sequence) GetCharge() *int {
	return s.charge
}

// SetCharge sets the charge state
func (s *Sequence) SetCharge(charge *int) {
	s.charge = charge
}

// GetIonicSpecies returns the ionic species
func (s *Sequence) GetIonicSpecies() *string {
	return s.ionicSpecies
}

// SetIonicSpecies sets the ionic species
func (s *Sequence) SetIonicSpecies(species *string) {
	s.ionicSpecies = species
}

// IsChimeric returns whether the sequence is chimeric
func (s *Sequence) IsChimeric() bool {
	return s.isChimeric
}

// GetPeptidoforms returns the peptidoforms for chimeric sequences
func (s *Sequence) GetPeptidoforms() []*Sequence {
	return s.peptidoforms
}

// IsMultiChain returns whether the sequence is multi-chain
func (s *Sequence) IsMultiChain() bool {
	return s.isMultiChain
}

// GetChains returns the chains for multi-chain sequences
func (s *Sequence) GetChains() []*Sequence {
	return s.chains
}