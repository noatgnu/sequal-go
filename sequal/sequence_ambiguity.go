package sequal

import "fmt"

// SequenceAmbiguity represents ambiguity in the amino acid sequence
type SequenceAmbiguity struct {
	Value    string
	Position int
}

// NewSequenceAmbiguity creates a new SequenceAmbiguity instance
func NewSequenceAmbiguity(value string, position int) *SequenceAmbiguity {
	return &SequenceAmbiguity{
		Value:    value,
		Position: position,
	}
}

// GetValue returns the ambiguous sequence value
func (sa *SequenceAmbiguity) GetValue() string {
	return sa.Value
}

// GetPosition returns the position of the ambiguity
func (sa *SequenceAmbiguity) GetPosition() int {
	return sa.Position
}

// String returns a string representation of the sequence ambiguity
func (sa *SequenceAmbiguity) String() string {
	return fmt.Sprintf("SequenceAmbiguity(value='%s', position=%d)", sa.Value, sa.Position)
}
