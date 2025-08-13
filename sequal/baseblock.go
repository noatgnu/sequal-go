package sequal

import (
	"crypto/sha256"
	"encoding/json"
	"fmt"
)

// BaseBlock represents a biochemical building block with position and mass properties.
// This is the base interface for components like amino acids, modifications, etc.
type BaseBlock interface {
	GetValue() string
	GetPosition() *int
	SetPosition(position *int)
	IsBranch() bool
	GetMass() *float64
	SetMass(mass *float64)
	GetExtra() interface{}
	SetExtra(value interface{})
	ToMap() map[string]interface{}
	String() string
	Equal(other BaseBlock) bool
	Hash() (string, error)
}

// BaseBlockImpl provides a concrete implementation of the BaseBlock interface
type BaseBlockImpl struct {
	value    string
	position *int
	branch   bool
	mass     *float64
	extra    interface{}
}

// NewBaseBlock creates a new BaseBlock instance
func NewBaseBlock(value string, position *int, branch bool, mass *float64) *BaseBlockImpl {
	return &BaseBlockImpl{
		value:    value,
		position: position,
		branch:   branch,
		mass:     mass,
		extra:    nil,
	}
}

// GetValue returns the identifier of the block
func (b *BaseBlockImpl) GetValue() string {
	return b.value
}

// GetPosition returns the position of the block
func (b *BaseBlockImpl) GetPosition() *int {
	return b.position
}

// SetPosition sets the position of the block
func (b *BaseBlockImpl) SetPosition(position *int) {
	b.position = position
}

// IsBranch checks if the block is a branch
func (b *BaseBlockImpl) IsBranch() bool {
	return b.branch
}

// GetMass returns the mass of the block
func (b *BaseBlockImpl) GetMass() *float64 {
	return b.mass
}

// SetMass sets the mass of the block
func (b *BaseBlockImpl) SetMass(mass *float64) {
	b.mass = mass
}

// GetExtra returns extra information associated with the block
func (b *BaseBlockImpl) GetExtra() interface{} {
	return b.extra
}

// SetExtra sets extra information for the block
func (b *BaseBlockImpl) SetExtra(value interface{}) {
	b.extra = value
}

// ToMap converts the block to a map representation
func (b *BaseBlockImpl) ToMap() map[string]interface{} {
	return map[string]interface{}{
		"value":    b.value,
		"position": b.position,
		"branch":   b.branch,
		"mass":     b.mass,
		"extra":    b.extra,
	}
}

// String returns a string representation of the block
func (b *BaseBlockImpl) String() string {
	return b.value
}

// Equal checks if two blocks are equal
func (b *BaseBlockImpl) Equal(other BaseBlock) bool {
	if other == nil {
		return false
	}
	bHash, err := b.Hash()
	if err != nil {
		return false
	}
	otherHash, err := other.Hash()
	if err != nil {
		return false
	}
	return bHash == otherHash
}

func (b *BaseBlockImpl) Hash() (string, error) {
	/// calculate has of the block based on ToMap
	var bMap = b.ToMap()
	jsonData, err := json.Marshal(bMap)
	if err != nil {
		return "", err
	}
	hash := sha256.Sum256(jsonData)
	return fmt.Sprintf("%x", hash), nil
}

// Helper function to compare pointers
func equalPtr(a, b *int) bool {
	if a == nil && b == nil {
		return true
	}
	if a == nil || b == nil {
		return false
	}
	return *a == *b
}
