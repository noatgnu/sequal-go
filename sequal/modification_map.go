package sequal

// ModificationMap maps modifications to their positions in a sequence for quick lookup
type ModificationMap struct {
	seq             string
	ignorePositions map[int]bool
	modDictByName   map[string]*Modification
	modPositionDict map[string][]int
	positionToMods  map[int][]*Modification
}

// NewModificationMap creates a new ModificationMap instance
func NewModificationMap(
	seq string,
	mods []*Modification,
	ignorePositions map[int]bool,
	parsePosition bool,
	modPositionDict map[string][]int,
) *ModificationMap {
	if ignorePositions == nil {
		ignorePositions = make(map[int]bool)
	}

	if modPositionDict == nil {
		modPositionDict = make(map[string][]int)
	}

	mm := &ModificationMap{
		seq:             seq,
		ignorePositions: ignorePositions,
		modDictByName:   make(map[string]*Modification),
		modPositionDict: modPositionDict,
		positionToMods:  make(map[int][]*Modification),
	}

	mm.buildMappings(mods, parsePosition)
	return mm
}

// buildMappings builds internal mappings between modifications and positions
func (mm *ModificationMap) buildMappings(mods []*Modification, parsePosition bool) {
	for _, mod := range mods {
		modName := mod.String()
		mm.modDictByName[modName] = mod

		if parsePosition {
			if mod.GetRegex() != nil && mod.GetPosition() == nil {
				for _, match := range mod.FindPositions(mm.seq) {
					startPos, endPos := match[0], match[1]
					for pos := startPos; pos < endPos; pos++ {
						if !mm.ignorePositions[pos] {
							if _, exists := mm.positionToMods[pos]; !exists {
								mm.positionToMods[pos] = []*Modification{}
							}
							mm.positionToMods[pos] = append(mm.positionToMods[pos], mod)

							if _, exists := mm.modPositionDict[modName]; !exists {
								mm.modPositionDict[modName] = []int{}
							}
							mm.modPositionDict[modName] = append(mm.modPositionDict[modName], pos)
						}
					}
				}
			} else if mod.GetPosition() != nil {
				pos := *mod.GetPosition()
				if !mm.ignorePositions[pos] {
					if _, exists := mm.positionToMods[pos]; !exists {
						mm.positionToMods[pos] = []*Modification{}
					}
					mm.positionToMods[pos] = append(mm.positionToMods[pos], mod)

					if _, exists := mm.modPositionDict[modName]; !exists {
						mm.modPositionDict[modName] = []int{}
					}
					mm.modPositionDict[modName] = append(mm.modPositionDict[modName], pos)
				}
			}
		}
	}
}

// GetModPositions gets the positions of a modification by its name
func (mm *ModificationMap) GetModPositions(modName string) []int {
	positions, exists := mm.modPositionDict[modName]
	if !exists {
		return nil
	}
	return positions
}

// GetMod gets the Modification object by its name
func (mm *ModificationMap) GetMod(modName string) *Modification {
	mod, exists := mm.modDictByName[modName]
	if !exists {
		return nil
	}
	return mod
}

// GetModsAtPosition gets all modifications at a specific position
func (mm *ModificationMap) GetModsAtPosition(position int) []*Modification {
	mods, exists := mm.positionToMods[position]
	if !exists {
		return []*Modification{}
	}
	return mods
}

// HasModAtPosition checks if a position has any modification or a specific modification
func (mm *ModificationMap) HasModAtPosition(position int, modName *string) bool {
	mods := mm.GetModsAtPosition(position)
	if len(mods) == 0 {
		return false
	}

	if modName == nil {
		return len(mods) > 0
	}

	for _, mod := range mods {
		if mod.String() == *modName {
			return true
		}
	}
	return false
}

// ToDict converts the modification map to a dictionary representation
func (mm *ModificationMap) ToDict() map[string]interface{} {
	modifications := make(map[string]interface{})

	for modName, positions := range mm.modPositionDict {
		mod := mm.GetMod(modName)
		modifications[modName] = map[string]interface{}{
			"positions": positions,
			"details":   mod.ToMap(),
		}
	}

	return map[string]interface{}{
		"sequence":      mm.seq,
		"modifications": modifications,
	}
}
