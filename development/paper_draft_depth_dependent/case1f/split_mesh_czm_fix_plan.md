# Fix: Pre-split Mesh + CZM Interface Compatibility

## Problem Statement

When using `--split-mesh` / `--use-split` with meshes that contain CZM (Cohesive Zone
Model) interfaces created by `BreakMeshByBlockGenerator`, the solve fails with:

```
Element 60125 on side 0 is missing a neighbor (hence identified as an external side)
but has interface kernel(s) defined on the boundary.
```

## Root Cause

**The disjoint neighbor boundary pair information is NOT serialized by libMesh's
CheckpointIO.**

### Data Flow

1. **`BreakMeshByBlockGenerator`** duplicates nodes along block interfaces to create
   CZM boundaries. It then calls:
   ```cpp
   mesh.add_disjoint_neighbor_boundary_pairs(boundary_id, reverse_boundary_id,
                                              RealVectorValue(0, 0, 0));
   ```
   This registers pairs of boundary IDs as "disjoint neighbors" in
   `MeshBase::_disjoint_neighbor_boundary_pairs` (a `PeriodicBoundaries` object).
   libMesh uses this to know that elements across these boundaries should be ghosted
   to each other even though they no longer share nodes.

2. **`SplitMeshAction`** calls `libMesh::split_mesh()` then writes via `CheckpointIO`.
   The mesh is fully prepared at this point (disjoint pairs are set, ghosting works).

3. **`CheckpointIO::write()`** serializes: nodes, elements, subdomain names, boundary
   IDs/names, nodesets, extra integers. It does **NOT** serialize
   `_disjoint_neighbor_boundary_pairs`. (Confirmed: no reference to
   `disjoint_neighbor` anywhere in `checkpoint_io.C` or its header.)

4. **`--use-split` read phase**: `SetupMeshAction` (line 314) skips the entire
   MeshGenerator system, so `BreakMeshByBlockGenerator` never runs and its
   relationship manager (`ElementSideNeighborLayers`) is never registered.
   `CheckpointIO::read()` reconstructs the mesh geometry but without the disjoint
   neighbor pairs.

5. **Assembly phase**: `InterfaceKernelBase` expects neighbor elements across the CZM
   boundary. Since ghosting wasn't set up (no disjoint neighbor info), the neighbor
   element is on another processor and not ghosted. Error is thrown in
   `NonlinearThread.C` line 285.

## Proposed Fixes

Two complementary fixes at different levels:

### Fix 1: libMesh CheckpointIO (root cause fix)

Serialize `_disjoint_neighbor_boundary_pairs` in CheckpointIO so the data survives
write/read cycles.

**Files to modify:**
- `libmesh/src/mesh/checkpoint_io.C`
- `libmesh/include/libmesh/checkpoint_io.h`

#### In `CheckpointIO::write()` — add after boundary info writing:

```cpp
// Write disjoint neighbor boundary pairs (needed for CZM interfaces)
#ifdef LIBMESH_ENABLE_PERIODIC
{
  const auto * disjoint_pairs = mesh.get_disjoint_neighbor_boundary_pairs();
  std::vector<boundary_id_type> pair_data;  // [b1, b2, b1, b2, ...]
  std::vector<Real> translation_data;       // [tx, ty, tz, tx, ty, tz, ...]

  if (disjoint_pairs)
  {
    // PeriodicBoundaries is a std::map<boundary_id_type, PeriodicBoundaryBase*>
    for (const auto & [bid, pb] : *disjoint_pairs)
    {
      pair_data.push_back(bid);
      pair_data.push_back(pb->pairedboundary);
      const auto & t = pb->get_transformation_matrix();
      // For CZM, translation is stored; extract from the periodic boundary
      // For simplicity, store the translation vector
      translation_data.push_back(pb->get_translation_vector()(0));
      translation_data.push_back(pb->get_translation_vector()(1));
      translation_data.push_back(pb->get_translation_vector()(2));
    }
  }

  // Write count + data
  io.data(pair_data, "# disjoint neighbor boundary pairs");
  io.data(translation_data, "# disjoint neighbor translations");
}
#endif
```

#### In `CheckpointIO::read()` — add after boundary info reading:

```cpp
// Read disjoint neighbor boundary pairs
#ifdef LIBMESH_ENABLE_PERIODIC
{
  std::vector<boundary_id_type> pair_data;
  std::vector<Real> translation_data;

  if (io.data_exists("# disjoint neighbor boundary pairs"))
  {
    io.data(pair_data, "# disjoint neighbor boundary pairs");
    io.data(translation_data, "# disjoint neighbor translations");

    for (std::size_t i = 0; i < pair_data.size(); i += 2)
    {
      RealVectorValue trans(translation_data[3*(i/2)],
                            translation_data[3*(i/2)+1],
                            translation_data[3*(i/2)+2]);
      mesh.add_disjoint_neighbor_boundary_pairs(pair_data[i], pair_data[i+1], trans);
    }
  }
}
#endif
```

**Note:** The exact serialization API depends on whether CheckpointIO uses Xdr or
another format. The above is pseudocode — the actual implementation needs to follow
CheckpointIO's existing patterns for `io.data()` calls.

### Fix 2: MOOSE-level workaround (if libMesh fix is too invasive)

Store the disjoint neighbor boundary pair info in MOOSE's mesh metadata (which IS
written during split), and restore it when reading split files.

**Files to modify:**
- `framework/src/meshgenerators/BreakMeshByBlockGenerator.C`
- `framework/src/actions/SetupMeshAction.C` (or `MooseMesh.C`)

#### Step A: Save pairs as mesh metadata in BreakMeshByBlockGenerator

In `BreakMeshByBlockGenerator::generate()`, after creating all interface boundaries,
store the boundary pair info as mesh metadata:

```cpp
// In BreakMeshByBlockGenerator::generate(), after the addInterfaceBoundary loop:

// Save disjoint neighbor pairs as mesh metadata for split mesh compatibility
std::vector<boundary_id_type> disjoint_pair_ids;
for (const auto & [subid_pair, bid] : _subid_pairs_to_boundary_id)
{
  auto rev_pair = std::make_pair(subid_pair.second, subid_pair.first);
  if (_subid_pairs_to_boundary_id.count(rev_pair))
  {
    disjoint_pair_ids.push_back(bid);
    disjoint_pair_ids.push_back(_subid_pairs_to_boundary_id[rev_pair]);
  }
}
declareMeshProperty("disjoint_neighbor_boundary_pairs", disjoint_pair_ids);
```

#### Step B: Restore pairs when reading split mesh

In `SetupMeshAction::act()` for the `init_mesh` task, after the mesh is read from
split files, check for the metadata and re-apply:

```cpp
// In SetupMeshAction::act(), init_mesh task, after _mesh->init():

if (_use_split)
{
  // Restore disjoint neighbor boundary pairs from mesh metadata
  // (needed for CZM interfaces created by BreakMeshByBlockGenerator)
  const auto & mg_names = _app.getMeshGeneratorNames();
  for (const auto & mg_name : mg_names)
  {
    if (hasMeshProperty<std::vector<boundary_id_type>>(
            "disjoint_neighbor_boundary_pairs", mg_name))
    {
      const auto & pairs = getMeshProperty<std::vector<boundary_id_type>>(
          "disjoint_neighbor_boundary_pairs", mg_name);
      auto & lm_mesh = _mesh->getMesh();
      for (std::size_t i = 0; i + 1 < pairs.size(); i += 2)
        lm_mesh.add_disjoint_neighbor_boundary_pairs(
            pairs[i], pairs[i + 1], RealVectorValue(0, 0, 0));
    }
  }
}
```

## Recommended Approach

**Fix 2 (MOOSE-level)** is recommended as the first step because:
1. It doesn't require a libMesh PR and rebuild
2. MOOSE's mesh metadata IS already written during split (`SplitMeshAction` line 97-101
   writes `writeRestartableMetaData(MESH_META_DATA, ...)`)
3. Mesh metadata IS restored when reading split files (`FileMesh.C` line 153:
   `possiblyLoadRestartableMetaData(MESH_META_DATA, ...)`)
4. It's self-contained within MOOSE framework code

**Fix 1 (libMesh-level)** is the long-term robust solution and should be submitted as a
libMesh PR separately.

## Key Source Files Reference

| File | Role |
|------|------|
| `BreakMeshByBlockGenerator.C:602` | Calls `add_disjoint_neighbor_boundary_pairs()` |
| `SplitMeshAction.C:90` | Calls `libMesh::split_mesh()` |
| `SplitMeshAction.C:97-101` | Writes mesh metadata (survives split!) |
| `SetupMeshAction.C:314` | Skips MeshGenerators when `_use_split` |
| `FileMesh.C:153` | Loads mesh metadata when reading checkpoint |
| `checkpoint_io.C` | Does NOT serialize disjoint pairs |
| `InterfaceKernelBase.C:68-72` | Declares ghosting RM |
| `NonlinearThread.C:285` | Throws "missing neighbor" error |
| `MeshBase.h:1996-2071` | Stores `_disjoint_neighbor_boundary_pairs` |
