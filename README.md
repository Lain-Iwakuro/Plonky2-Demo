## Understand Plonky2

### Example: Matrix Multiplication

我们通过一个矩阵乘法的样例电路来理解 `Plonky2` 的工作方式，源文件位于`plonky2/src/bin`，`main.rs` 中的代码如下:

```rust
use anyhow::Result;
use plonky2::field::types::Field;
use plonky2::iop::target::Target;
use plonky2::iop::witness::{PartialWitness, WitnessWrite};
use plonky2::plonk::circuit_builder::CircuitBuilder;
use plonky2::plonk::circuit_data::CircuitConfig;
use plonky2::plonk::config::{GenericConfig, PoseidonGoldilocksConfig};
// Randomness generation.
use rand_chacha::ChaChaRng;
use rand_chacha::rand_core::SeedableRng;
use rand::Rng;

use rayon::ThreadPoolBuilder;

/// An example of using Plonky2 to prove a statement of the form
/// "I know A * B = C".
fn main() -> Result<()> {
    // Use only 2 threads.
    ThreadPoolBuilder::new().num_threads(1).build_global().unwrap();

    const D: usize = 2;
    type C = PoseidonGoldilocksConfig;
    type F = <C as GenericConfig<D>>::F;

    let config = CircuitConfig::standard_recursion_config();
    let mut builder = CircuitBuilder::<F, D>::new(config);

    // The arithmetic circuit.
    // let m = 64;
    let m = 2;
    let mut a = Vec::<Vec<Target>>::new();
    let mut b = Vec::<Vec<Target>>::new();
    let mut c = Vec::<Vec<Target>>::new();
    for i in 0..m {
        a.push(Vec::<Target>::new());
        b.push(Vec::<Target>::new());
        c.push(Vec::<Target>::new());
        (0..m).for_each(|_j: usize| {
            a[i].push(builder.add_virtual_target());
            b[i].push(builder.add_virtual_target());
        });
    }
    for i in 0..m {
        for j in 0..m {
            let initial = builder.constant(F::from_canonical_u32(0));
            let mut current_target = initial;
            for k in 0..m {
                let c_i_j_k = builder.mul(a[i][k], b[k][j]);
                current_target = builder.add(current_target, c_i_j_k);
            }
            c[i].push(current_target);
        }
    }

    // Public inputs are the initial value (provided below) and the result (which is generated).
    for i in 0..m {
        for j in 0..m {
            builder.register_public_input(a[i][j]);
            builder.register_public_input(b[i][j]);
            builder.register_public_input(c[i][j]);
        }
    }
    
    // Build the circuit without actually calculating the result.
    let data = builder.build::<C>();
    
    // let mut rand = ChaChaRng::from_entropy();
    let mut pw = PartialWitness::new();
    let inputs: Vec<Vec<Vec<u32>>> = vec![vec![vec![1, 2], vec![3, 4]], vec![vec![5, 6], vec![7, 8]]];
    for i in 0..m {
        for j in 0..m {
            // pw.set_target(a[i][j], F::from_canonical_u32(rand.gen_range(u32::MIN..u32::MAX)));
            pw.set_target(a[i][j], F::from_canonical_u32(inputs[0][i][j]));
            // pw.set_target(b[i][j], F::from_canonical_u32(rand.gen_range(u32::MIN..u32::MAX)));
            pw.set_target(b[i][j], F::from_canonical_u32(inputs[1][i][j]));
        }
    }
    
    // Construct the proof.
    let proof = data.prove(pw)?;

    // Show the length of the public inputs of proof (should equal 3 * (m ** 2)).
    println!("length of proof.public_inputs is {}", proof.public_inputs.len());
    
    // Verify the proof.
    data.verify(proof)
}

```



### Arithmetic

在这一部分中，我们创建一个 `CircuitBuilder` 结构体的实例 `builder`，再调用它的函数为电路添加组件. 

#### Target, Virtual Target and Wire

首先为 A B 两个矩阵创建 `Target`，这里调用了 `builder.add_virtual_target`，它的作用是为 `builder` 添加一个 `virtual_target `:

```rust
let config = CircuitConfig::standard_recursion_config();
let mut builder = CircuitBuilder::<F, D>::new(config);

// The arithmetic circuit.
// let m = 64;
let m = 2;
let mut a = Vec::<Vec<Target>>::new();
let mut b = Vec::<Vec<Target>>::new();
let mut c = Vec::<Vec<Target>>::new();
for i in 0..m {
    a.push(Vec::<Target>::new());
    b.push(Vec::<Target>::new());
    c.push(Vec::<Target>::new());
    (0..m).for_each(|_j: usize| {
        a[i].push(builder.add_virtual_target());
        b[i].push(builder.add_virtual_target());
    });
}
```
```rust
/// Adds a new "virtual" target. This is not an actual wire in the witness, but just a target
/// that help facilitate witness generation. In particular, a generator can assign a values to a
/// virtual target, which can then be copied to other (virtual or concrete) targets. When we
/// generate the final witness (a grid of wire values), these virtual targets will go away.
pub fn add_virtual_target(&mut self) -> Target {
    let index = self.virtual_target_index;
    self.virtual_target_index += 1;
    Target::VirtualTarget { index }
}
```

`VirtualTarget` 是一种 `Target`，它在电路中并没有实际的位置，只包含它自己的 `index`，每添加一个 `index` 就加一，从 0 开始；另一种 `Target` 是 `Wire`，包含自己的行 `row` 和列 `column` :

```rust
pub enum Target {
    Wire(Wire),
    /// A target that doesn't have any inherent location in the witness (but it can be copied to
    /// another target that does). This is useful for representing intermediate values in witness
    /// generation.
    VirtualTarget {
        index: usize,
    },
}

pub struct Wire {
    /// Row index of the wire.
    pub row: usize,
    /// Column index of the wire.
    pub column: usize,
}
```

**注意**: `VirtualTarget` 和 `Wire` 都不包含任何实际数值，只是对电路不同位置的编号.

#### Gates

接着，为了计算矩阵乘法，我们需要计算许多乘累加，这通过调用 `builder.mul` 和 `builder.add` 实现:

```rust
for i in 0..m {
    for j in 0..m {
        let initial = builder.constant(F::from_canonical_u32(0));
        let mut current_target = initial;
        for k in 0..m {
            let c_i_j_k = builder.mul(a[i][k], b[k][j]);
            current_target = builder.add(current_target, c_i_j_k);
        }
        c[i].push(current_target);
    }
}
```

##### Constants

不过我们首先要将累加的初始值设为 0，通过调用 `builder.constant` 实现. 这个函数会先判断常数 `c` 是否已经添加过，如果没有添加过则调用 `self.add_virtual_target` :

```rust
/// Returns a routable target with the given constant value.
pub fn constant(&mut self, c: F) -> Target {
    if let Some(&target) = self.constants_to_targets.get(&c) {
        // We already have a wire for this constant.
        return target;
    }

    let target = self.add_virtual_target();
    self.constants_to_targets.insert(c, target);
    self.targets_to_constants.insert(target, c);

    target
}
```

##### Multiply and Add

之后对于每次运算，调用 `builder.mul` 或 `builder.add`，它们都会调用 `self.arithmetic` 再接着调用 `self.add_base_arithmetic_operation`，因为它们都可以表达为 `BaseArithmeticOperation` 即 `result = const_0 * multiplicand_0 * multiplicand_1 + const_1 * addend` 的形式，乘法门 `const_0 = 1, const_1 = 1`，加法门 `const_0 = const_1 = 1, multiplicand_1 = 1 ` :

```rust
/// Computes `x * y`.
pub fn mul(&mut self, x: Target, y: Target) -> Target {
    // x * y = 1 * x * y + 0 * x
    self.arithmetic(F::ONE, F::ZERO, x, y, x)
}
/// Computes `x + y`.
pub fn add(&mut self, x: Target, y: Target) -> Target {
    let one = self.one();
    // x + y = 1 * x * 1 + 1 * y
    self.arithmetic(F::ONE, F::ONE, x, one, y)
}

pub(crate) struct BaseArithmeticOperation<F: Field64> {
    const_0: F,
    const_1: F,
    multiplicand_0: Target,
    multiplicand_1: Target,
    addend: Target,
}
```

在 `add_base_arithmetic_operation` 中，调用 `self.find_slot` 函数为新的运算寻找位置，得到行数 `gate` 和 该运算的序号 `i`，各个运算数的位置分别为 `(gate, 4i), (gate, 4i+1), (gate, 4i+2)`，运算结果的位置为 `(gate, 4i+3)`. 还需要调用 `self.connect` 把取值相同位置不同的 `target` 绑定起来.

```rust
fn add_base_arithmetic_operation(&mut self, operation: BaseArithmeticOperation<F>) -> Target {
    let gate = ArithmeticGate::new_from_config(&self.config);
    let constants = vec![operation.const_0, operation.const_1];
    let (gate, i) = self.find_slot(gate, &constants, &constants);
    let wires_multiplicand_0 = Target::wire(gate, ArithmeticGate::wire_ith_multiplicand_0(i));
    let wires_multiplicand_1 = Target::wire(gate, ArithmeticGate::wire_ith_multiplicand_1(i));
    let wires_addend = Target::wire(gate, ArithmeticGate::wire_ith_addend(i));

    self.connect(operation.multiplicand_0, wires_multiplicand_0);
    self.connect(operation.multiplicand_1, wires_multiplicand_1);
    self.connect(operation.addend, wires_addend);

    Target::wire(gate, ArithmeticGate::wire_ith_output(i))
}
```

##### Find Slot

一行同一类型的门（如 `const_0 = 1, const_1 = 0` 的乘法门）能容纳 `num_ops` 个运算 ，在 `find_slot` 中，先判断是否有与待添加门同一类型且未填满的 `slot`，如果有，则将 `(gate_idx, slot_idx) = slot` 作为返回值用以添加运算；否则，调用 `self.add_gate` 添加新的一行（一组新的同类型的门）.

当一整个 `slot` 填满时，把相应 `gate` 的参数从 `current_slot` 中除去，这样下次出现相同运算时仍然需要重新调用 `self.add_gate`.

```rust
/// Find an available slot, of the form `(row, op)` for gate `G` using parameters `params`
/// and constants `constants`. Parameters are any data used to differentiate which gate should be
/// used for the given operation.
pub fn find_slot<G: Gate<F, D> + Clone>(
    &mut self,
    gate: G,
    params: &[F],
    constants: &[F],
) -> (usize, usize) {
    let num_gates = self.num_gates();
    let num_ops = gate.num_ops();
    let gate_ref = GateRef::new(gate.clone());
    let gate_slot = self.current_slots.entry(gate_ref.clone()).or_default();
    let slot = gate_slot.current_slot.get(params);
    let (gate_idx, slot_idx) = if let Some(&s) = slot {
        s
    } else {
        self.add_gate(gate, constants.to_vec());
        (num_gates, 0)
    };
    let current_slot = &mut self.current_slots.get_mut(&gate_ref).unwrap().current_slot;
    if slot_idx == num_ops - 1 {
        // We've filled up the slots at this index.
        current_slot.remove(params);
    } else {
        // Increment the slot operation index.
        current_slot.insert(params.to_vec(), (gate_idx, slot_idx + 1));
    }

    (gate_idx, slot_idx)
}
```

`add_gate` 函数会检查 `gate_type` 是否合法，扩充 `self.constant_generators` （算术门不需要），并将 `gate_ref` 和 `constants` 作为 `GateInstance` 加进 `self.gate_instances` 里，并返回行号.

`GateRef` 具体的运作方式==有待研究==.

```rust
/// Adds a gate to the circuit, and returns its index.
pub fn add_gate<G: Gate<F, D>>(&mut self, gate_type: G, mut constants: Vec<F>) -> usize {
    self.check_gate_compatibility(&gate_type);

    assert!(
        constants.len() <= gate_type.num_constants(),
        "Too many constants."
    );
    constants.resize(gate_type.num_constants(), F::ZERO);

    let row = self.gate_instances.len();

    self.constant_generators
    .extend(gate_type.extra_constant_wires().into_iter().map(
        |(constant_index, wire_index)| ConstantGenerator {
            row,
            constant_index,
            wire_index,
            constant: F::ZERO, // Placeholder; will be replaced later.
        },
    ));

    // Note that we can't immediately add this gate's generators, because the list of constants
    // could be modified later, i.e. in the case of `ConstantGate`. We will add them later in
    // `build` instead.

    // Register this gate type if we haven't seen it before.
    let gate_ref = GateRef::new(gate_type);
    self.gates.insert(gate_ref.clone());

    self.gate_instances.push(GateInstance {
        gate_ref,
        constants,
    });

    row
}
```

##### Connect

`connect` 会将两个取值相等且 `is_routable` 的 `target` 绑定加入 `self.copy_constraints`，在后续计算 `wire_partition` 和生成 `witness` 时起作用.

```rust
/// Uses Plonk's permutation argument to require that two elements be equal.
/// Both elements must be routable, otherwise this method will panic.
pub fn connect(&mut self, x: Target, y: Target) {
    assert!(
        x.is_routable(&self.config),
        "Tried to route a wire that isn't routable"
    );
    assert!(
        y.is_routable(&self.config),
        "Tried to route a wire that isn't routable"
    );
    self.copy_constraints
    .push(CopyConstraint::new((x, y), self.context_log.open_stack()));
}
```

#### Public Inputs

调用 `builder.register_public_input` 来把电路输入和计算结果设为 Public Inputs.

```rust
// Public inputs are the initial value (provided below) and the result (which is generated).
for i in 0..m {
    for j in 0..m {
        builder.register_public_input(a[i][j]);
        builder.register_public_input(b[i][j]);
        builder.register_public_input(c[i][j]);
    }
}

/// Registers the given target as a public input.
pub fn register_public_input(&mut self, target: Target) {
    self.public_inputs.push(target);
}
```



### Build

在这一部分中，main 调用 `CircuitBuilder` 结构体的函数  `build()` 得到一个 `CircuitData` 结构体的实例 `data`.

`build` 函数的声明:

```rust
pub fn build<C: GenericConfig<D, F = F>>(mut self) -> CircuitData<F, C, D>;
```

而 `CircuitData` 结构体包含以下成员:

```rust
pub struct CircuitData<F: RichField + Extendable<D>, C: GenericConfig<D, F = F>, const D: usize> {
    pub prover_only: ProverOnlyCircuitData<F, C, D>,
    pub verifier_only: VerifierOnlyCircuitData<C, D>,
    pub common: CommonCircuitData<F, D>,
}
```

其中的 `verifier_only` 包含 `constants_sigmas_cap` 和 `circuit_digest`，可以理解为对电路结构的承诺.

接着来看 `build` 函数都做了什么来生成 `CircuitData`.

#### More Gates

- 为电路的 Public Inputs 添加 `PoseidonGate`
- 添加 `PublicInputGate` 并将 `PoseidonGate` 的输出与 `PublicInputGate` 绑定
- 为每一个不使用的 `PublicInputGate` 的端口添加 `RandomValueGenerator`
- 添加 `ConstantGate`，将其的端口与常数对象绑定，为排序后的每一个常数创建并添加 `ConstantGenerator`

```rust
// Hash the public inputs, and route them to a `PublicInputGate` which will enforce that
// those hash wires match the claimed public inputs.
let num_public_inputs = self.public_inputs.len();
let public_inputs_hash =
    self.hash_n_to_hash_no_pad::<C::InnerHasher>(self.public_inputs.clone()); // 添加 PoseidonGate
let pi_gate = self.add_gate(PublicInputGate, vec![]);
for (&hash_part, wire) in public_inputs_hash
    .elements
    .iter()
    .zip(PublicInputGate::wires_public_inputs_hash())
{
    self.connect(hash_part, Target::wire(pi_gate, wire)) // 将 public_inputs_hash 与 PublicInputGate 绑定
}
self.randomize_unused_pi_wires(pi_gate); // 添加 RandomValueGenerator
// Constant Gates and Generators...
```

##### Poseidon Gate

对 Public Inputs 按照 8 个一组（最后一组可能少于 8 个）分组，每一组对应 1 个 `PoseidonGate`.

一个 `PoseidonGate` 有 12 个输入端口 (0~11) 和 12 个输出端口 (12~23)，以及一个 `swap_wire` 端口 (24) 与常数 0 绑定. 其中前 8 个（或不足 8 个）输入端口与电路的 Public Input 绑定，剩余输入端口与前一个 `PoseidonGate` 的相应输出端口绑定（例如，(3, 8) ~ (3, 11) 依次与 (2, 20) ~ (2, 23) 绑定）. 对于第一个 `PoseidonGate`，剩余输入端口与 `target` 常数 0 绑定.

`self.hash_n_to_hash_no_pad` 调用 `hash_n_to_m_no_pad` :（最后一段的 Squeeze 步骤是为了之后与 `PublicInputGate` 连接）

```rust
pub fn hash_n_to_m_no_pad<H: AlgebraicHasher<F>>(
    &mut self,
    inputs: Vec<Target>,
    num_outputs: usize,
) -> Vec<Target> {
    let zero = self.zero();
    let mut state = H::AlgebraicPermutation::new(std::iter::repeat(zero)); // 初始状态为 0

    // Absorb all input chunks.
    for input_chunk in inputs.chunks(H::AlgebraicPermutation::RATE) {
        // Overwrite the first r elements with the inputs. This differs from a standard sponge,
        // where we would xor or add in the inputs. This is a well-known variant, though,
        // sometimes called "overwrite mode".
        state.set_from_slice(input_chunk, 0); // 用 inputs 覆盖新 PoseidonGate 的一部分输入端口
        state = self.permute::<H>(state); // 添加 PoseidonGate，其输出端口作为新的输入端口
    }

    // Squeeze until we have the desired number of outputs.
    let mut outputs = Vec::with_capacity(num_outputs);
    loop {
        for &s in state.squeeze() {
            outputs.push(s);
            if outputs.len() == num_outputs { // 本例中 num_outputs = 4，不需要 loop
                return outputs;
            }
        }
        state = self.permute::<H>(state);
    }
}
```

`self.permute` 调用 `permute_swapped ` :

```rust
fn permute_swapped<const D: usize>(
    inputs: Self::AlgebraicPermutation,
    swap: BoolTarget,
    builder: &mut CircuitBuilder<F, D>,
) -> Self::AlgebraicPermutation
where
F: RichField + Extendable<D>,
{
    let gate_type = PoseidonGate::<F, D>::new();
    let gate = builder.add_gate(gate_type, vec![]);

    let swap_wire = PoseidonGate::<F, D>::WIRE_SWAP;
    let swap_wire = Target::wire(gate, swap_wire);
    builder.connect(swap.target, swap_wire);

    // Route input wires.
    let inputs = inputs.as_ref();
    for i in 0..SPONGE_WIDTH {
        let in_wire = PoseidonGate::<F, D>::wire_input(i);
        let in_wire = Target::wire(gate, in_wire);
        builder.connect(inputs[i], in_wire); // 绑定输入端口
    }

    // Collect output wires.
    Self::AlgebraicPermutation::new(
        (0..SPONGE_WIDTH).map(|i| Target::wire(gate, PoseidonGate::<F, D>::wire_output(i))),
    )
}
```

##### Public Input Gate

通过 `self.add_gate` 添加，绑定了最后一个 `PoseidonGate` 的输出. 它似乎没有计算上的作用，其他效果==有待研究==.

对于剩余端口，调用 `self.randomize_unused_pi_wires`，为 `PublicInputGate` 的每个剩余端口添加一个 `RandomValueGenerator`.

```rust
/// In PLONK's permutation argument, there's a slight chance of division by zero. We can
/// mitigate this by randomizing some unused witness elements, so if proving fails with
/// division by zero, the next attempt will have an (almost) independent chance of success.
/// See https://github.com/mir-protocol/plonky2/issues/456
fn randomize_unused_pi_wires(&mut self, pi_gate: usize) {
    for wire in PublicInputGate::wires_public_inputs_hash().end..self.config.num_wires {
        self.add_simple_generator(RandomValueGenerator {
            target: Target::wire(pi_gate, wire),
        });
    }
}
```

##### Constant Gate

以下代码先添加了 `ConstantGate`，再为排序后的每个常数生成并添加了一个 `ConstantGenerator`.

```rust
// Make sure we have enough constant generators. If not, add a `ConstantGate`.
while self.constants_to_targets.len() > self.constant_generators.len() {
    self.add_gate(
        ConstantGate {
            num_consts: self.config.num_constants,
        },
        vec![],
    );
}

// For each constant-target pair used in the circuit, use a constant generator to fill this target.
for ((c, t), mut const_gen) in self
    .constants_to_targets
    .clone()
    .into_iter()
    // We need to enumerate constants_to_targets in some deterministic order to ensure that
    // building a circuit is deterministic.
    .sorted_by_key(|(c, _t)| c.to_canonical_u64())
    .zip(self.constant_generators.clone())
{
    // Set the constant in the constant polynomial.
    self.gate_instances[const_gen.row].constants[const_gen.constant_index] = c;
    // Generate a copy between the target and the routable wire.
    self.connect(Target::wire(const_gen.row, const_gen.wire_index), t);
    // Set the constant in the generator (it's initially set with a dummy value).
    const_gen.set_constant(c);
    self.add_simple_generator(const_gen);
}
```

#### Add Generators

- 为实际需要运算的门的每一次运算添加一个 generator
- 生成表达计算依赖关系的 `generator_indices_by_watches`

##### Incomplete Gates

每个 `ArithmeticGate` 占据一整行，可以容纳 `num_ops` 个同类运算，但行中可能有空位，被称为 incomplete gates，下面的代码获取了每行实际结束的位置.

```rust
// Map between gates where not all generators are used and the gate's number of used generators.
let incomplete_gates = self
    .current_slots
    .values()
    .flat_map(|current_slot| current_slot.current_slot.values().copied())
    .collect::<HashMap<_, _>>();
```

为每次运算生成 generator 时，去除多余的 generators.

```rust
self.add_generators(
    self.gate_instances
    .iter()
    .enumerate()
    .flat_map(|(index, gate)| {
        let mut gens = gate.gate_ref.0.generators(index, &gate.constants);
        // Remove unused generators, if any.
        if let Some(&op) = incomplete_gates.get(&index) {
            gens.drain(op..);
        }
        gens
    })
    .collect(),
);
```

调用 `gate.gate_ref.0.generators` 得到（可能很多个） generators，`generators` 函数会按 gate 类型不同返回不同结果.

然后调用 `self.add_generators` 添加 generators:

```rust
pub fn add_generators(&mut self, generators: Vec<WitnessGeneratorRef<F>>) {
    self.generators.extend(generators);
}
```

##### Watch Lists

实际电路的运算存在数据的依赖关系，所以把每个 generator（即每次计算）所需要的 `Vec<Target>` 记录下来:

```rust
// Index generator indices by their watched targets.
let mut generator_indices_by_watches = BTreeMap::new();
for (i, generator) in self.generators.iter().enumerate() {
    for watch in generator.0.watch_list() {
        let watch_index = forest.target_index(watch);
        let watch_rep_index = forest.parents[watch_index];
        generator_indices_by_watches
        .entry(watch_rep_index)
        .or_insert_with(Vec::new)
        .push(i);
    }
}
```

（此处用到 `forest`，是在计算完电路划分后执行的）

其中调用 `generator.0.watch_list` 来获取依赖的 `Vec<Target>`.

```rust
/// Targets to be "watched" by this generator. Whenever a target in the watch list is populated,
/// the generator will be queued to run.
fn watch_list(&self) -> Vec<Target> {
    self.inner.dependencies()
}
fn dependencies(&self) -> Vec<Target>;
```

不同的 generator 有不同的 `dependencies` 实现，如 `ArithmeticBaseGenerator` 和 `PoseidonGenerator` 的实现分别为:

```rust
fn dependencies(&self) -> Vec<Target> {
    [
        ArithmeticGate::wire_ith_multiplicand_0(self.i),
        ArithmeticGate::wire_ith_multiplicand_1(self.i),
        ArithmeticGate::wire_ith_addend(self.i),
    ]
    .iter()
    .map(|&i| Target::wire(self.row, i))
    .collect()
}
fn dependencies(&self) -> Vec<Target> {
    (0..SPONGE_WIDTH)
    .map(|i| PoseidonGate::<F, D>::wire_input(i))
    .chain(Some(PoseidonGate::<F, D>::WIRE_SWAP))
    .map(|column| Target::wire(self.row, column))
    .collect()
}
```

#### Pre-compute Things

- 计算常数向量和选择器多项式
- 计算电路划分和 sigma 向量
- 计算 FFT table
- 生成有关电路结构的承诺

##### Constant Vector and Selector Polynomial

Plonk 协议需要用多项式来表示门的属性，即 `selector_polynomials`. 而每个多项式的度数不能超过 `max_degree`，因此把每一类门按`gate.0.degree`排序，用贪心算法分组，然后得到一个反映分组信息的 `selector_indices` 和一个 `groups.len() * max_degree` 的 `polynomials: Vec<PolynomialValues<F>>`.

==需要再研究 `size` 和 `gate.0.degree` 为什么要加起来和 `max_degree` 比较==.

```rust
let quotient_degree_factor = self.config.max_quotient_degree_factor;
let mut gates = self.gates.iter().cloned().collect::<Vec<_>>();
// Gates need to be sorted by their degrees (and ID to make the ordering deterministic) to compute the selector polynomials.
gates.sort_unstable_by_key(|g| (g.0.degree(), g.0.id()));
let (mut constant_vecs, selectors_info) =
	selector_polynomials(&gates, &self.gate_instances, quotient_degree_factor + 1);

/// Returns the selector polynomials and related information.
///
/// Selector polynomials are computed as follows:
/// Partition the gates into (the smallest amount of) groups `{ G_i }`, such that for each group `G`
/// `|G| + max_{g in G} g.degree() <= max_degree`. These groups are constructed greedily from
/// the list of gates sorted by degree.
/// We build a selector polynomial `S_i` for each group `G_i`, with
/// S_i\[j\] =
///     if j-th row gate=g_k in G_i
///         k
///     else
///         UNUSED_SELECTOR
pub(crate) fn selector_polynomials<F: RichField + Extendable<D>, const D: usize>(
    gates: &[GateRef<F, D>],
    instances: &[GateInstance<F, D>],
    max_degree: usize,
) -> (Vec<PolynomialValues<F>>, SelectorsInfo);
```
##### Subgroup, Coset Shifts and Sigma Vector

Plonk 协议计算 sigma 向量需要一个 `F` 上阶为 `degree` 的乘法子群，调用`two_adic_subgroup` 和 `primitive_rout_of_unity` 来计算，原理是将最原始的单位根平方 `TWO_ADICITY - n_log` 次.

```rust
let subgroup = F::two_adic_subgroup(degree_bits);

/// Computes the subgroup generated by the root of unity of a given order generated by `Self::primitive_root_of_unity`.
fn two_adic_subgroup(n_log: usize) -> Vec<Self> {
    let generator = Self::primitive_root_of_unity(n_log);
    generator.powers().take(1 << n_log).collect()
}

fn primitive_root_of_unity(n_log: usize) -> Self {
    assert!(n_log <= Self::TWO_ADICITY);
    let base = Self::POWER_OF_TWO_GENERATOR;
    base.exp_power_of_2(Self::TWO_ADICITY - n_log)
}
```
此外还要计算子群的 `self.config.num_routed_wires` 个陪集移位，取生成元的 0~num 次方即可.

```rust
let k_is = get_unique_coset_shifts(degree, self.config.num_routed_wires);

/// Finds a set of shifts that result in unique cosets for the multiplicative subgroup of size
/// `2^subgroup_bits`.
pub fn get_unique_coset_shifts<F: Field>(subgroup_size: usize, num_shifts: usize) -> Vec<F> {
    // From Lagrange's theorem.
    let num_cosets = (F::order() - 1u32) / (subgroup_size as u32);
    assert!(
        BigUint::from(num_shifts) <= num_cosets,
        "The subgroup does not have enough distinct cosets"
    );

    // Let g be a generator of the entire multiplicative group. Let n be the order of the subgroup.
    // The subgroup can be written as <g^(|F*| / n)>. We can use g^0, ..., g^(num_shifts - 1) as our
    // shifts, since g^i <g^(|F*| / n)> are distinct cosets provided i < |F*| / n, which we checked.
    F::MULTIPLICATIVE_GROUP_GENERATOR // i.e. 7
        .powers()
        .take(num_shifts)
        .collect()
}
```

之后调用 `self.sigma_vec` 计算 `sigma_vecs` 和电路划分，其中的 `forest` 是一个并查集，利用 `self.copy_constraints` 来合并绑定了的 target，调用 `forest.compress_paths` 来把树的高度降为 1.

```rust
let (sigma_vecs, forest) = timed!(
    timing,
    "generate sigma polynomials",
    self.sigma_vecs(&k_is, &subgroup)
);

fn sigma_vecs(&self, k_is: &[F], subgroup: &[F]) -> (Vec<PolynomialValues<F>>, Forest) {
    let degree = self.gate_instances.len();
    let degree_log = log2_strict(degree);
    let config = &self.config;
    let mut forest = Forest::new(
        config.num_wires,
        config.num_routed_wires,
        degree,
        self.virtual_target_index,
    );

    for gate in 0..degree {
        for input in 0..config.num_wires {
            forest.add(Target::Wire(Wire {
                row: gate,
                column: input,
            }));
        }
    }

    for index in 0..self.virtual_target_index {
        forest.add(Target::VirtualTarget { index });
    }

    for &CopyConstraint { pair: (a, b), .. } in &self.copy_constraints {
        forest.merge(a, b);
    }

    forest.compress_paths();

    let wire_partition = forest.wire_partition();
    (
        wire_partition.get_sigma_polys(degree_log, k_is, subgroup),
        forest,
    )
}
```

调用 `forest.wire_partition` 来生成 `wire_partition`，即一个值为向量的哈希表，包含了每个根结点的全部叶子结点

```rust
/// Assumes `compress_paths` has already been called.
pub fn wire_partition(&mut self) -> WirePartition {
    let mut partition = HashMap::<_, Vec<_>>::new();

    // Here we keep just the Wire targets, filtering out everything else.
    for row in 0..self.degree {
        for column in 0..self.num_routed_wires {
            let w = Wire { row, column };
            let t = Target::Wire(w);
            let x_parent = self.parents[self.target_index(t)];
            partition.entry(x_parent).or_default().push(w);
        }
    }

    let partition = partition.into_values().collect();
    WirePartition { partition }
}
```

最后调用 `wire_partition.get_sigma_polys` 来计算 sigma 多项式的值:

先调用 `self.get_sigma_map` 得到 sigma 映射，再生成 `num_routed_wires` 个长为 `degree` 的 `PolynomialValues`，值为 `k_is[x / degree] * subgroup[x % degree]`.

因为 `subgroup` 的大小只为 `degree`，而实际 `wire_values` 大小还需乘上 `num_wires`，所以需要 `k_is` 来扩大集合的元素个数. 设 `subgroup` 生成元为 $h$，`k_is` 生成元为 $g$，则若位置为 $(i,j)$ 被映射到位置 $(i',j')$，则 `sigma_polys` 在 $(i,j)$ 处的取值就为 $g^{j'}h^{i'}$. 在 `prove` 函数中会用到这一点，参见 Prove-Openings 小节.

```rust
pub(crate) fn get_sigma_polys<F: Field>(
    &self,
    degree_log: usize,
    k_is: &[F],
    subgroup: &[F],
) -> Vec<PolynomialValues<F>> {
    let degree = 1 << degree_log;
    let sigma = self.get_sigma_map(degree, k_is.len());

    sigma
    .chunks(degree)
    .map(|chunk| {
        let values = chunk
        .par_iter()
        .map(|&x| k_is[x / degree] * subgroup[x % degree])
        .collect::<Vec<_>>();
        PolynomialValues::new(values)
    })
    .collect()
}
```

 sigma 映射将每个集合的 wire 循环映射到自身:

比如，若原始电路划分为 $\{\{a,b,c\},\{d,e\},\{f\}\}$，则 $[a,\dots,f]$ 的映射结果为 $[b,c,a,e,d,f]$.

```rust
/// Generates sigma in the context of Plonk, which is a map from `[kn]` to `[kn]`, where `k` is
/// the number of routed wires and `n` is the number of gates.
fn get_sigma_map(&self, degree: usize, num_routed_wires: usize) -> Vec<usize> {
    // Find a wire's "neighbor" in the context of Plonk's "extended copy constraints" check. In
    // other words, find the next wire in the given wire's partition. If the given wire is last in
    // its partition, this will loop around. If the given wire has a partition all to itself, it
    // is considered its own neighbor.
    let mut neighbors = HashMap::new();
    for subset in &self.partition {
        for n in 0..subset.len() {
            neighbors.insert(subset[n], subset[(n + 1) % subset.len()]);
        }
    }

    let mut sigma = Vec::new();
    for column in 0..num_routed_wires {
        for row in 0..degree {
            let wire = Wire { row, column };
            let neighbor = neighbors[&wire];
            sigma.push(neighbor.column * degree + neighbor.row);
        }
    }
    sigma
}
```

##### FFT Table

计算 FFT 运算所需要的根.计算`constants_sigmas_commitment` 时使用.

```rust
// Precompute FFT roots.
let max_fft_points = 1 << (degree_bits + max(rate_bits, log2_ceil(quotient_degree_factor)));
let fft_root_table = fft_root_table(max_fft_points);

pub fn fft_root_table<F: Field>(n: usize) -> FftRootTable<F> {
    let lg_n = log2_strict(n);
    // bases[i] = g^2^i, for i = 0, ..., lg_n - 1
    let mut bases = Vec::with_capacity(lg_n);
    let mut base = F::primitive_root_of_unity(lg_n);
    bases.push(base);
    for _ in 1..lg_n {
        base = base.square(); // base = g^2^_
        bases.push(base);
    }

    let mut root_table = Vec::with_capacity(lg_n);
    for lg_m in 1..=lg_n {
        let half_m = 1 << (lg_m - 1);
        let base = bases[lg_n - lg_m];
        let root_row = base.powers().take(half_m.max(2)).collect();
        root_table.push(root_row);
    }
    root_table
}
```

##### Commitment of Circuit

调用 `PolynomialBatch::from_values` 计算 `constants_sigmas_commitment` :

```rust
let constants_sigmas_vecs = [constant_vecs, sigma_vecs.clone()].concat();
let constants_sigmas_commitment = PolynomialBatch::<F, C, D>::from_values(
    constants_sigmas_vecs,
    rate_bits,
    PlonkOracle::CONSTANTS_SIGMAS.blinding,
    cap_height,
    &mut timing,
    Some(&fft_root_table),
);
```

先 IFFT 再对系数调用 `from_coeffs` :

```rust
/// Creates a list polynomial commitment for the polynomials interpolating the values in `values`.
pub fn from_values(
    values: Vec<PolynomialValues<F>>,
    rate_bits: usize,
    blinding: bool,
    cap_height: usize,
    timing: &mut TimingTree,
    fft_root_table: Option<&FftRootTable<F>>,
) -> Self {
    let coeffs = timed!(
        timing,
        "IFFT",
        values.into_par_iter().map(|v| v.ifft()).collect::<Vec<_>>()
    );

    Self::from_coeffs(
        coeffs,
        rate_bits,
        blinding,
        cap_height,
        timing,
        fft_root_table,
    )
}
```

`from_coeffs` 调用 `lde_values` 来计算 low degree extension 的值，再调用 `MerkleTree::new` 来建树

```rust
/// Creates a list polynomial commitment for the polynomials `polynomials`.
pub fn from_coeffs(
    polynomials: Vec<PolynomialCoeffs<F>>,
    rate_bits: usize,
    blinding: bool,
    cap_height: usize,
    timing: &mut TimingTree,
    fft_root_table: Option<&FftRootTable<F>>,
) -> Self {
    let degree = polynomials[0].len();
    let lde_values = timed!(
        timing,
        "FFT + blinding",
        Self::lde_values(&polynomials, rate_bits, blinding, fft_root_table)
    );

    let mut leaves = timed!(timing, "transpose LDEs", transpose(&lde_values));
    reverse_index_bits_in_place(&mut leaves);
    let merkle_tree = timed!(
        timing,
        "build Merkle tree",
        MerkleTree::new(leaves, cap_height)
    );

    Self {
        polynomials,
        merkle_tree,
        degree_log: log2_strict(degree),
        rate_bits,
        blinding,
    }
}
```

`lde_values` 调用 `PolynomialCoeffs::lde` 将系数向量零填充至长度为原来的 `2 ** rate_bits` 倍

```rust
fn lde_values(
    polynomials: &[PolynomialCoeffs<F>],
    rate_bits: usize,
    blinding: bool,
    fft_root_table: Option<&FftRootTable<F>>,
) -> Vec<Vec<F>> {
    let degree = polynomials[0].len();

    // If blinding, salt with two random elements to each leaf vector.
    let salt_size = if blinding { SALT_SIZE } else { 0 };

    polynomials
    .par_iter()
    .map(|p| {
        assert_eq!(p.len(), degree, "Polynomial degrees inconsistent");
        p.lde(rate_bits)
        .coset_fft_with_options(F::coset_shift(), Some(rate_bits), fft_root_table)
        .values
    })
    .chain(
        (0..salt_size)
        .into_par_iter()
        .map(|_| F::rand_vec(degree << rate_bits)),
    )
    .collect()
}
```

再调用 `coset_fft_with_options` 计算 FFT，得到多项式扩展后的值.

```rust
/// Returns the evaluation of the polynomial on the coset `shift*H`.
pub fn coset_fft_with_options(
    &self,
    shift: F,
    zero_factor: Option<usize>,
    root_table: Option<&FftRootTable<F>>,
) -> PolynomialValues<F> {
    let modified_poly: Self = shift
    .powers()
    .zip(&self.coeffs)
    .map(|(r, &c)| r * c)
    .collect::<Vec<_>>()
    .into();
    modified_poly.fft_with_options(zero_factor, root_table)
}

pub fn fft_with_options(
    self,
    zero_factor: Option<usize>,
    root_table: Option<&FftRootTable<F>>,
) -> PolynomialValues<F> {
    fft_with_options(self, zero_factor, root_table)
}

#[inline]
pub fn fft_with_options<F: Field>(
    poly: PolynomialCoeffs<F>,
    zero_factor: Option<usize>,
    root_table: Option<&FftRootTable<F>>,
) -> PolynomialValues<F> {
    let PolynomialCoeffs { coeffs: mut buffer } = poly;
    fft_dispatch(&mut buffer, zero_factor, root_table); // 这里面求了 FFT
    PolynomialValues::new(buffer)
}
```

`circuit_digest` 体现了电路总体的结构信息，通过调用 `hash_no_pad` 求几部分的哈希得到.

```rust
let constants_sigmas_cap = constants_sigmas_commitment.merkle_tree.cap.clone();
let domain_separator = self.domain_separator.unwrap_or_default();
let domain_separator_digest = C::Hasher::hash_pad(&domain_separator);
let circuit_digest_parts = [
    constants_sigmas_cap.flatten(),
    domain_separator_digest.to_vec(),
    vec![
        F::from_canonical_usize(degree_bits),
        /* Add other circuit data here */
    ],
];
let circuit_digest = C::Hasher::hash_no_pad(&circuit_digest_parts.concat());
```



### Prove

在这一部分，我们设置两个矩阵的值分别为 $\begin{bmatrix}1 & 2 \\ 3 & 4\end{bmatrix}$ 和 $\begin{bmatrix}5 & 6 \\ 7 & 8\end{bmatrix}$，计算 `witness` 并调用 `data.prove` 生成 `proof`.

```rust
let mut pw = PartialWitness::new();
let inputs: Vec<Vec<Vec<u32>>> = vec![vec![vec![1, 2], vec![3, 4]], vec![vec![5, 6], vec![7, 8]]];
for i in 0..m {
    for j in 0..m {
        pw.set_target(a[i][j], F::from_canonical_u32(inputs[0][i][j]));
        pw.set_target(b[i][j], F::from_canonical_u32(inputs[1][i][j]));
    }
}

// Construct the proof.
let proof = data.prove(pw)?;
```

`proof` 为 `ProofWithPublicInputs` 结构体的实例，包含一个 `proof` 和一个 `public_inputs`.

```rust
pub struct ProofWithPublicInputs<
    F: RichField + Extendable<D>,
    C: GenericConfig<D, F = F>,
    const D: usize,
> {
    pub proof: Proof<F, C, D>,
    pub public_inputs: Vec<F>,
}
```

而 `Proof` 结构体有以下成员:

```rust
pub struct Proof<F: RichField + Extendable<D>, C: GenericConfig<D, F = F>, const D: usize> {
    /// Merkle cap of LDEs of wire values.
    pub wires_cap: MerkleCap<F, C::Hasher>,
    /// Merkle cap of LDEs of Z, in the context of Plonk's permutation argument.
    pub plonk_zs_partial_products_cap: MerkleCap<F, C::Hasher>,
    /// Merkle cap of LDEs of the quotient polynomial components.
    pub quotient_polys_cap: MerkleCap<F, C::Hasher>,
    /// Purported values of each polynomial at the challenge point.
    pub openings: OpeningSet<F, D>,
    /// A batch FRI argument for all openings.
    pub opening_proof: FriProof<F, C::Hasher, D>,
}
```

接着来看 `ProofWithPublicInputs` 是如何生成的.

#### Witness Generation

##### Partial Witness

`PartialWitness` 包含一个 `HashMap<Target, F>`，在调用 `prove` 函数之前，需要先调用 `pw.set_target` 为电路的输入赋值，每个 `Target` 至多赋一个值:

```rust
let mut pw = PartialWitness::new();
let inputs: Vec<Vec<Vec<u32>>> = vec![vec![vec![1, 2], vec![3, 4]], vec![vec![5, 6], vec![7, 8]]];
for i in 0..m {
    for j in 0..m {
        pw.set_target(a[i][j], F::from_canonical_u32(inputs[0][i][j]));
        pw.set_target(b[i][j], F::from_canonical_u32(inputs[1][i][j]));
    }    
}

pub struct PartialWitness<F: Field> {
    pub(crate) target_values: HashMap<Target, F>,
}

fn set_target(&mut self, target: Target, value: F) {
    let opt_old_value = self.target_values.insert(target, value);
    if let Some(old_value) = opt_old_value {
        assert_eq!(
            value, old_value,
            "Target {:?} was set twice with different values: {} != {}",
            target, old_value, value
        );
    }
}
```

##### Partition Witness

然后在 `main.rs` 中调用 `data.prove(pw)`，继而调用同名函数:

```rust
pub fn prove(&self, inputs: PartialWitness<F>) -> Result<ProofWithPublicInputs<F, C, D>> {
    prove::<F, C, D>(
        &self.prover_only,
        &self.common,
        inputs,
        &mut TimingTree::default(),
    )
}

pub fn prove<F: RichField + Extendable<D>, C: GenericConfig<D, F = F>, const D: usize>(
    prover_data: &ProverOnlyCircuitData<F, C, D>,
    common_data: &CommonCircuitData<F, D>,
    inputs: PartialWitness<F>,
    timing: &mut TimingTree,
) -> Result<ProofWithPublicInputs<F, C, D>>
where
    C::Hasher: Hasher<F>,
    C::InnerHasher: Hasher<F>;
```

在 `prove` 中，首先调用 `generate_partial_witness` 生成 `partition_witness`，`PartitionWitness` 结构体包含一个并查集每个集合的取值 `values`，每个集合中元素（`target`）索引到代表性元素索引的映射 `representative_mat` 以及一些电路信息:

```rust
let mut partition_witness = timed!(
    timing,
    &format!("run {} generators", prover_data.generators.len()),
    generate_partial_witness(inputs, prover_data, common_data)
);

/// `PartitionWitness` holds a disjoint-set forest of the targets respecting a circuit's copy constraints.
/// The value of a target is defined to be the value of its root in the forest.
#[derive(Clone)]
pub struct PartitionWitness<'a, F: Field> {
    pub values: Vec<Option<F>>,
    pub representative_map: &'a [usize],
    pub num_wires: usize,
    pub degree: usize,
}
```

`generate_partial_witness` 首先调用 `witness.set_target` 将部分 `target` 设置 `inputs: PartialWitness` 即 `pw` 中提供的电路输入，再生成一个包含所有 `generators` 索引的列表 `pending_generator_indices`（初始长度为 `generators` 个数），之后对 `pending_generator_indices` 中所有 `generator` 依次调用 `generator.0.run`，如果运算成功执行则将这个 `generator` 标记；

将运算输出在 `buffer.target_values` 里的值取出，通过调用 `witness.set_target_returning_rep` 将其加入 `witness` 并得到 `target` 在 `wire_partition` 中的代表性元素，之后基于 `build` 时得到的 `generator_indices_by_watch` 将那些依赖本次运算结果且未成功执行过的 `generator` 的索引加入 `pending_generator_indices`，让它们在下一轮循环中尝试运算，以保证最终所有 `generator` 得以成功执行.

不过，如果 `generators` 中的 `generator` 本来就是按计算依赖顺序出现的，==似乎==并不存在无法成功执行的情况，此时循环只需要一轮就能成功执行所有运算，此处代码的设计似乎只是为了保险.

```rust
/// Given a `PartitionWitness` that has only inputs set, populates the rest of the witness using the
/// given set of generators.
pub(crate) fn generate_partial_witness<
    'a,
    F: RichField + Extendable<D>,
    C: GenericConfig<D, F = F>,
    const D: usize,
>(
    inputs: PartialWitness<F>,
    prover_data: &'a ProverOnlyCircuitData<F, C, D>,
    common_data: &'a CommonCircuitData<F, D>,
) -> PartitionWitness<'a, F> {
    let config = &common_data.config;
    let generators = &prover_data.generators;
    let generator_indices_by_watches = &prover_data.generator_indices_by_watches;

    let mut witness = PartitionWitness::new(
        config.num_wires,
        common_data.degree(),
        &prover_data.representative_map,
    );

    for (t, v) in inputs.target_values.into_iter() {
        witness.set_target(t, v);
    }

    // Build a list of "pending" generators which are queued to be run. Initially, all generators
    // are queued.
    let mut pending_generator_indices: Vec<_> = (0..generators.len()).collect();

    // We also track a list of "expired" generators which have already returned false.
    let mut generator_is_expired = vec![false; generators.len()];
    let mut remaining_generators = generators.len();

    let mut buffer = GeneratedValues::empty();

    // Keep running generators until we fail to make progress.
    while !pending_generator_indices.is_empty() {
        let mut next_pending_generator_indices = Vec::new();

        for &generator_idx in &pending_generator_indices {
            if generator_is_expired[generator_idx] {
                continue;
            }

            let finished = generators[generator_idx].0.run(&witness, &mut buffer);
            if finished {
                generator_is_expired[generator_idx] = true;
                remaining_generators -= 1;
            }

            // Merge any generated values into our witness, and get a list of newly-populated
            // targets' representatives.
            let new_target_reps = buffer
                .target_values
                .drain(..)
                .flat_map(|(t, v)| witness.set_target_returning_rep(t, v));

            // Enqueue unfinished generators that were watching one of the newly populated targets.
            for watch in new_target_reps {
                let opt_watchers = generator_indices_by_watches.get(&watch);
                if let Some(watchers) = opt_watchers {
                    for &watching_generator_idx in watchers {
                        if !generator_is_expired[watching_generator_idx] {
                            next_pending_generator_indices.push(watching_generator_idx);
                        }
                    }
                }
            }
        }

        pending_generator_indices = next_pending_generator_indices;
    }

    assert_eq!(
        remaining_generators, 0,
        "{} generators weren't run",
        remaining_generators,
    );

    witness
}
```

`set_target` 也调用了 `set_target_returning_rep`，把值赋给 `target` 所在的 `partition` 的代表性元素并返回其索引.

```rust
fn set_target(&mut self, target: Target, value: F) {
    self.set_target_returning_rep(target, value);
}

/// Set a `Target`. On success, returns the representative index of the newly-set target. If the
/// target was already set, returns `None`.
pub(crate) fn set_target_returning_rep(&mut self, target: Target, value: F) -> Option<usize> {
    let rep_index = self.representative_map[self.target_index(target)];
    let rep_value = &mut self.values[rep_index];
    if let Some(old_value) = *rep_value {
        assert_eq!(
            value, old_value,
            "Partition containing {:?} was set twice with different values: {} != {}",
            target, old_value, value
        );
        None
    } else {
        *rep_value = Some(value);
        Some(rep_index)
    }
}
```

`generator.0.run` 对于不同类型 `generator` 的实现也不同，但它们都会将计算结果放入 `out_buffer` 并返回计算是否成功.

```rust
/// Run this generator, returning a flag indicating whether the generator is finished. If the
/// flag is true, the generator will never be run again, otherwise it will be queued for another
/// run next time a target in its watch list is populated.
fn run(&self, witness: &PartitionWitness<F>, out_buffer: &mut GeneratedValues<F>) -> bool {
    if witness.contains_all(&self.inner.dependencies()) {
        self.inner.run_once(witness, out_buffer);
        true
    } else {
        false
    }
}

fn run_once(&self, witness: &PartitionWitness<F>, out_buffer: &mut GeneratedValues<F>);

// ArithmeticBaseGenerator
fn run_once(&self, witness: &PartitionWitness<F>, out_buffer: &mut GeneratedValues<F>) {
    let get_wire = |wire: usize| -> F { witness.get_target(Target::wire(self.row, wire)) };

    let multiplicand_0 = get_wire(ArithmeticGate::wire_ith_multiplicand_0(self.i));
    let multiplicand_1 = get_wire(ArithmeticGate::wire_ith_multiplicand_1(self.i));
    let addend = get_wire(ArithmeticGate::wire_ith_addend(self.i));

    let output_target = Target::wire(self.row, ArithmeticGate::wire_ith_output(self.i));

    let computed_output =
    multiplicand_0 * multiplicand_1 * self.const_0 + addend * self.const_1;

    out_buffer.set_target(output_target, computed_output)
}

// RandomValueGenerator
fn run_once(&self, _witness: &PartitionWitness<F>, out_buffer: &mut GeneratedValues<F>) {
    let random_value = F::rand();
    out_buffer.set_target(self.target, random_value);
}
```

##### Full Witness

`PartialWitness` 只包含电路输入值，而 `PartitionWitness` 中每个 `partition` 只有一个 `target` 被赋值，为了得到每个 `target` 都被赋值的 witness，调用 `full_witness` 得到 `num_wires * degree` 大小的 `MatrixWitness`，其每一列对于一个门:

```rust
let witness = timed!(
    timing,
    "compute full witness",
    partition_witness.full_witness()
);

pub fn full_witness(self) -> MatrixWitness<F> {
    let mut wire_values = vec![vec![F::ZERO; self.degree]; self.num_wires];
    for i in 0..self.degree {
        for j in 0..self.num_wires {
            let t = Target::Wire(Wire { row: i, column: j });
            if let Some(x) = self.try_get_target(t) {
                wire_values[j][i] = x;
            }
        }
    }

    MatrixWitness { wire_values }
}
```

##### Commitment of Witness

`MatrixWitness` 很大，需要转变为 `wires_commitment` 以便后续处理，这里同样调用了`PolynomialBatch::from_values`:

```rust
let wires_commitment = timed!(
    timing,
    "compute wires commitment",
    PolynomialBatch::<F, C, D>::from_values(
        wires_values,
        config.fri_config.rate_bits,
        config.zero_knowledge && PlonkOracle::WIRES.blinding,
        config.fri_config.cap_height,
        timing,
        prover_data.fft_root_table.as_ref(),
    )
);
```

#### Challenges

得到了 witness 之后，证明者需要接受来自验证者的随机挑战，但来回通信并不方便，这里使用了 Fiat-Shamir 变换将挑战转为由证明者自己生成，需要证明者计算包含电路结构、电路输入、计算结构的哈希，并将输出处理为来自验证者的挑战，由于哈希函数的性质，这样的挑战是近似随机的，也即不受证明者控制.

求 `public_inputs_hash` 仍是通过调用 `hash_no_pad` 实现的，之后调用 `challenger.observe_hash` 和 `challenger.observe_cap` 将信息交给哈希函数，再调用 `challenger.get_n_challenges` 获取一系列伪随机数:

```rust
let num_challenges = config.num_challenges;

let public_inputs = partition_witness.get_targets(&prover_data.public_inputs);
let public_inputs_hash = C::InnerHasher::hash_no_pad(&public_inputs);

let mut challenger = Challenger::<F, C::Hasher>::new();

// Observe the instance.
challenger.observe_hash::<C::Hasher>(prover_data.circuit_digest);
challenger.observe_hash::<C::InnerHasher>(public_inputs_hash);

challenger.observe_cap::<C::Hasher>(&wires_commitment.merkle_tree.cap);

let betas = challenger.get_n_challenges(num_challenges);
let gammas = challenger.get_n_challenges(num_challenges);
```

##### Observe

`observe_hash` 和 `observe_cap` 调用 `observe_elements`，将输入加入 `input_buffer` 每当 `input_buffer` 达到 `RATE` 时调用 `duplexing`，把 `input_buffer` 清空并把 `output_buffer` 填满:

```rust
pub fn observe_cap<OH: Hasher<F>>(&mut self, cap: &MerkleCap<F, OH>) {
    for &hash in &cap.0 {
        self.observe_hash::<OH>(hash);
    }
}

pub fn observe_hash<OH: Hasher<F>>(&mut self, hash: OH::Hash) {
    self.observe_elements(&hash.to_vec())
}

pub fn observe_elements(&mut self, elements: &[F]) {
    for &element in elements {
        self.observe_element(element);
    }
}

pub fn observe_element(&mut self, element: F) {
    // Any buffered outputs are now invalid, since they wouldn't reflect this input.
    self.output_buffer.clear();

    self.input_buffer.push(element);

    if self.input_buffer.len() == H::Permutation::RATE {
        self.duplexing();
    }
}
```

`duplexing` 将 `self.sponge_state` 的前几项用 `input_buffer` 替换，然后调用 `permute` 进而调用 `poseidon` 计算哈希并输出至 `output_buffer` :

```rust
/// Absorb any buffered inputs. After calling this, the input buffer will be empty, and the
/// output buffer will be full.
fn duplexing(&mut self) {
    assert!(self.input_buffer.len() <= H::Permutation::RATE);

    // Overwrite the first r elements with the inputs. This differs from a standard sponge,
    // where we would xor or add in the inputs. This is a well-known variant, though,
    // sometimes called "overwrite mode".
    self.sponge_state
    .set_from_iter(self.input_buffer.drain(..), 0);

    // Apply the permutation.
    self.sponge_state.permute();

    self.output_buffer.clear();
    self.output_buffer
    .extend_from_slice(self.sponge_state.squeeze());
}

fn permute(&mut self) {
    self.state = T::permute(self.state);
}

fn permute(input: [Self; SPONGE_WIDTH]) -> [Self; SPONGE_WIDTH] {
    <F as Poseidon>::poseidon(input)
}

fn poseidon(input: [Self; SPONGE_WIDTH]) -> [Self; SPONGE_WIDTH] {
    let mut state = input;
    let mut round_ctr = 0;

    Self::full_rounds(&mut state, &mut round_ctr);
    Self::partial_rounds(&mut state, &mut round_ctr);
    Self::full_rounds(&mut state, &mut round_ctr);
    debug_assert_eq!(round_ctr, N_ROUNDS);

    state
}
```

##### Get

`get_n_challenges` 调用 `get_challenge`，将 `challenger.output_buffer` 中的数取出，若`output_buffer` 为空则再调用 `duplexing` :

```rust
pub fn get_n_challenges(&mut self, n: usize) -> Vec<F> {
    (0..n).map(|_| self.get_challenge()).collect()
}

pub fn get_challenge(&mut self) -> F {
    // If we have buffered inputs, we must perform a duplexing so that the challenge will
    // reflect them. Or if we've run out of outputs, we must perform a duplexing to get more.
    if !self.input_buffer.is_empty() || self.output_buffer.is_empty() {
        self.duplexing();
    }

    self.output_buffer
    .pop()
    .expect("Output buffer should be non-empty")
}
```



#### Polynomials

现在需要按照 Plonk 协议的要求计算一些特定的多项式



##### Partial Products and Zs

调用 `all_wires_permutation_partial_products` 计算 `permutation polynomial` $Z(x)$ 和 `partial product polynomial` $\pi(x)$.

返回结果的最后 `num_challenges` 个向量为 $\beta,\gamma$ 不同时 $Z(x)$ 在 $H$ 上的取值，前面的为 $\pi(x)$，需要将其顺序调整一下，将其维度变为 `(num_challenges * num_routed_wires / 8, max_degree)`.

最后调用 `PolynomialBatch::from_values` 得到关于它们的 commitment `parial_products_zs_and_lookup_commitment` :

```rust
let mut partial_products_and_zs = timed!(
    timing,
    "compute partial products",
    all_wires_permutation_partial_products(&witness, &betas, &gammas, prover_data, common_data)
);

// Z is expected at the front of our batch; see `zs_range` and `partial_products_range`.
let plonk_z_vecs = partial_products_and_zs
    .iter_mut()
    .map(|partial_products_and_z| partial_products_and_z.pop().unwrap())
    .collect();
let zs_partial_products = [plonk_z_vecs, partial_products_and_zs.concat()].concat();

// All lookup polys: RE and partial SLDCs.
let lookup_polys =
compute_all_lookup_polys(&witness, &deltas, prover_data, common_data, has_lookup);

let zs_partial_products_lookups = if has_lookup {
    [zs_partial_products, lookup_polys].concat()
} else {
    zs_partial_products
};

let partial_products_zs_and_lookup_commitment = timed!(
    timing,
    "commit to partial products, Z's and, if any, lookup polynomials",
    PolynomialBatch::from_values(
        zs_partial_products_lookups,
        config.fri_config.rate_bits,
        config.zero_knowledge && PlonkOracle::ZS_PARTIAL_PRODUCTS.blinding,
        config.fri_config.cap_height,
        timing,
        prover_data.fft_root_table.as_ref(),
    )
);
```

若计算正确，有关系 $Z(x)\overset{r}{\underset{i=1}\Pi}{f_i'(x)}=Z(gx)\overset{r}{\underset{i=1}\Pi}{g_i'(x)},\forall{x\in{H=\{h,h^2,\dots,1\}}}$，且 $Z(h)=1$. 其中 $r$ 为 `num_routed_wires`，$f_i'(h^j)=f_i(h^j)+\beta\cdot{h^ig^j}+\gamma$，$g_i'(h^j)=f_i(h^j)+\beta\cdot{h^{i'}g^{j'}}+\gamma$，其中 $\sigma(i,j)=(i',j')$，$f_i(h^j)$ 为 $(i,j)$ 位置的取值.

由于 $r$ 很大，为降低约束的度数， Plonky2 将递推关系式按连乘分解为如下约束（为得到多项式形式需要把 $g'(x)$ 乘到左边）：

$\pi_1(x)=Z(x)\overset{8}{\underset{i=1}\Pi}{f_i'(x)/g_i'(x)},\pi_2(x)=\pi_1(x)\overset{16}{\underset{i=9}\Pi}{f_i'(x)/g_i'(x)},\dots,Z(gx)=\pi_s(x)\overset{r}{\underset{i=8s+1}\Pi}{f_i'(x)/g_i'(x)}$，其中 $s=\lfloor{r/8}\rfloor$.

在 `all_wires_permutation_partial_products` 中，对不同 `beta` 和 `gamma` 调用 `wires_permutation_partial_products_and_zs` 分开计算.

计算时先对不同 $i$，调用 `quotient_chunk_products` 8 个一组计算 $\Pi{f_i'(x)/g_i'(x)}$，递推地从 $h^0$ 到 $h^7$，从 $\pi_1$ 到 $\pi_9$ 到 $Z$ 计算取值，最后返回 10 个多项式在 $H$ 上的取值:

```rust
/// Compute the partial products used in the `Z` polynomials.
fn all_wires_permutation_partial_products<
    F: RichField + Extendable<D>,
    C: GenericConfig<D, F = F>,
    const D: usize,
>(
    witness: &MatrixWitness<F>,
    betas: &[F],
    gammas: &[F],
    prover_data: &ProverOnlyCircuitData<F, C, D>,
    common_data: &CommonCircuitData<F, D>,
) -> Vec<Vec<PolynomialValues<F>>> {
    (0..common_data.config.num_challenges)
        .map(|i| {
            wires_permutation_partial_products_and_zs(
                witness,
                betas[i],
                gammas[i],
                prover_data,
                common_data,
            )
        })
        .collect()
}

/// Compute the partial products used in the `Z` polynomial.
/// Returns the polynomials interpolating `partial_products(f / g)`
/// where `f, g` are the products in the definition of `Z`: `Z(g^i) = f / g`.
fn wires_permutation_partial_products_and_zs<
    F: RichField + Extendable<D>,
    C: GenericConfig<D, F = F>,
    const D: usize,
>(
    witness: &MatrixWitness<F>,
    beta: F,
    gamma: F,
    prover_data: &ProverOnlyCircuitData<F, C, D>,
    common_data: &CommonCircuitData<F, D>,
) -> Vec<PolynomialValues<F>> {
    let degree = common_data.quotient_degree_factor;
    let subgroup = &prover_data.subgroup;
    let k_is = &common_data.k_is;
    let num_prods = common_data.num_partial_products;
    let all_quotient_chunk_products = subgroup
        .par_iter()
        .enumerate()
        .map(|(i, &x)| {
            let s_sigmas = &prover_data.sigmas[i];
            let numerators = (0..common_data.config.num_routed_wires).map(|j| {
                let wire_value = witness.get_wire(i, j);
                let k_i = k_is[j];
                let s_id = k_i * x;
                wire_value + beta * s_id + gamma
            });
            let denominators = (0..common_data.config.num_routed_wires)
                .map(|j| {
                    let wire_value = witness.get_wire(i, j);
                    let s_sigma = s_sigmas[j];
                    wire_value + beta * s_sigma + gamma
                })
                .collect::<Vec<_>>();
            let denominator_invs = F::batch_multiplicative_inverse(&denominators);
            let quotient_values = numerators
                .zip(denominator_invs)
                .map(|(num, den_inv)| num * den_inv)
                .collect::<Vec<_>>();

            quotient_chunk_products(&quotient_values, degree)
        })
        .collect::<Vec<_>>();

    let mut z_x = F::ONE;
    let mut all_partial_products_and_zs = Vec::new();
    for quotient_chunk_products in all_quotient_chunk_products {
        let mut partial_products_and_z_gx =
            partial_products_and_z_gx(z_x, &quotient_chunk_products);
        // The last term is Z(gx), but we replace it with Z(x), otherwise Z would end up shifted.
        swap(&mut z_x, &mut partial_products_and_z_gx[num_prods]);
        all_partial_products_and_zs.push(partial_products_and_z_gx);
    }

    transpose(&all_partial_products_and_zs)
        .into_par_iter()
        .map(PolynomialValues::new)
        .collect()
}
```

`quotient_chunk_products` 按照 `max_degree` 计算连乘，将输入向量的长度缩小为 `max_degree` 分之一:

```rust
pub(crate) fn quotient_chunk_products<F: Field>(
    quotient_values: &[F],
    max_degree: usize,
) -> Vec<F> {
    debug_assert!(max_degree > 1);
    assert!(!quotient_values.is_empty());
    let chunk_size = max_degree;
    quotient_values
        .chunks(chunk_size)
        .map(|chunk| chunk.iter().copied().product())
        .collect()
}
```

`partial_products_and_z_gx` 计算向量 $[a_0,a_1,\dots,a_9]$ 的连乘，返回 $[a_0,a_0a_1,\dots,a_0a_1\cdots{a_9}]$ :

```rust
/// Compute partial products of the original vector `v` such that all products consist of `max_degree`
/// or less elements. This is done until we've computed the product `P` of all elements in the vector.
pub(crate) fn partial_products_and_z_gx<F: Field>(z_x: F, quotient_chunk_products: &[F]) -> Vec<F> {
    assert!(!quotient_chunk_products.is_empty());
    let mut res = Vec::new();
    let mut acc = z_x;
    for &quotient_chunk_product in quotient_chunk_products {
        acc *= quotient_chunk_product;
        res.push(acc);
    }
    res
}
```



##### Quotient Polynomials

之后再调用 `conpute_quotient_polys` 计算约束多项式 $C(x)$ 的系数，即 $C(x)$ 在低度扩展下的 `coset_fft`.

将其分为 `degree` 一组，再调用 `PolynomialBatch::from_coeffs` 来得到 commitment.

为什么这样分组==有待研究==.

```rust
challenger.observe_cap::<C::Hasher>(&partial_products_zs_and_lookup_commitment.merkle_tree.cap);

let alphas = challenger.get_n_challenges(num_challenges);

let quotient_polys = timed!(
    timing,
    "compute quotient polys",
    compute_quotient_polys::<F, C, D>(
        common_data,
        prover_data,
        &public_inputs_hash,
        &wires_commitment,
        &partial_products_zs_and_lookup_commitment,
        &betas,
        &gammas,
        &deltas,
        &alphas,
    )
);

let all_quotient_poly_chunks: Vec<PolynomialCoeffs<F>> = timed!(
    timing,
    "split up quotient polys",
    quotient_polys
    .into_par_iter()
    .flat_map(|mut quotient_poly| {
        quotient_poly.trim_to_len(quotient_degree).expect(
            "Quotient has failed, the vanishing polynomial is not divisible by Z_H",
        );
        // Split quotient into degree-n chunks.
        quotient_poly.chunks(degree)
    })
    .collect()
);

let quotient_polys_commitment = timed!(
    timing,
    "commit to quotient polys",
    PolynomialBatch::<F, C, D>::from_coeffs(
        all_quotient_poly_chunks,
        config.fri_config.rate_bits,
        config.zero_knowledge && PlonkOracle::QUOTIENT.blinding,
        config.fri_config.cap_height,
        timing,
        prover_data.fft_root_table.as_ref(),
    )
);
```

`compute_quotient_polys` 是一个复杂的函数，其主要思想是将约束 $c_i(x)$ 组装成 $C(x)=\Sigma_{i=0}^{l-1}\alpha^ic_i(x)$.

于是，为了计算 $C(x)$ 在各点的取值，需要先计算 $c_i(x)$. 在函数中，现将待求值的点（本例中为一个 64 阶乘法子群）分为 `BATCH_SIZE=32` 个一组，对每个 batch 提取各种信息，调用 `eval_vanishing_poly_base_batch` 求值，再除以 $g^8h^i-1,i=0\sim7$.

为什么这么除==有待研究==.

```rust
const BATCH_SIZE: usize = 32;

fn compute_quotient_polys<
    'a,
    F: RichField + Extendable<D>,
    C: GenericConfig<D, F = F>,
    const D: usize,
>(
    common_data: &CommonCircuitData<F, D>,
    prover_data: &'a ProverOnlyCircuitData<F, C, D>,
    public_inputs_hash: &<<C as GenericConfig<D>>::InnerHasher as Hasher<F>>::Hash,
    wires_commitment: &'a PolynomialBatch<F, C, D>,
    zs_partial_products_and_lookup_commitment: &'a PolynomialBatch<F, C, D>,
    betas: &[F],
    gammas: &[F],
    deltas: &[F],
    alphas: &[F],
) -> Vec<PolynomialCoeffs<F>> {
    let num_challenges = common_data.config.num_challenges;

    let has_lookup = common_data.num_lookup_polys != 0;

    let quotient_degree_bits = log2_ceil(common_data.quotient_degree_factor);
    assert!(
        quotient_degree_bits <= common_data.config.fri_config.rate_bits,
        "Having constraints of degree higher than the rate is not supported yet. \
        If we need this in the future, we can precompute the larger LDE before computing the `PolynomialBatch`s."
    );

    // We reuse the LDE computed in `PolynomialBatch` and extract every `step` points to get
    // an LDE matching `max_filtered_constraint_degree`.
    let step = 1 << (common_data.config.fri_config.rate_bits - quotient_degree_bits);
    // When opening the `Z`s polys at the "next" point in Plonk, need to look at the point `next_step`
    // steps away since we work on an LDE of degree `max_filtered_constraint_degree`.
    let next_step = 1 << quotient_degree_bits;

    let points = F::two_adic_subgroup(common_data.degree_bits() + quotient_degree_bits);
    let lde_size = points.len();

    let z_h_on_coset = ZeroPolyOnCoset::new(common_data.degree_bits(), quotient_degree_bits);

    let points_batches = points.par_chunks(BATCH_SIZE);
    let num_batches = ceil_div_usize(points.len(), BATCH_SIZE);
    let quotient_values: Vec<Vec<F>> = points_batches
        .enumerate()
        .flat_map(|(batch_i, xs_batch)| {
            // Each batch must be the same size, except the last one, which may be smaller.
            debug_assert!(
                xs_batch.len() == BATCH_SIZE
                    || (batch_i == num_batches - 1 && xs_batch.len() <= BATCH_SIZE)
            );

            let indices_batch: Vec<usize> =
                (BATCH_SIZE * batch_i..BATCH_SIZE * batch_i + xs_batch.len()).collect();

            let mut shifted_xs_batch = Vec::with_capacity(xs_batch.len());
            let mut local_zs_batch = Vec::with_capacity(xs_batch.len());
            let mut next_zs_batch = Vec::with_capacity(xs_batch.len());

            let mut local_lookup_batch = Vec::with_capacity(xs_batch.len());
            let mut next_lookup_batch = Vec::with_capacity(xs_batch.len());

            let mut partial_products_batch = Vec::with_capacity(xs_batch.len());
            let mut s_sigmas_batch = Vec::with_capacity(xs_batch.len());

            let mut local_constants_batch_refs = Vec::with_capacity(xs_batch.len());
            let mut local_wires_batch_refs = Vec::with_capacity(xs_batch.len());

            for (&i, &x) in indices_batch.iter().zip(xs_batch) {
                let shifted_x = F::coset_shift() * x;
                let i_next = (i + next_step) % lde_size;
                let local_constants_sigmas = prover_data
                    .constants_sigmas_commitment
                    .get_lde_values(i, step);
                let local_constants = &local_constants_sigmas[common_data.constants_range()];
                let s_sigmas = &local_constants_sigmas[common_data.sigmas_range()];
                let local_wires = wires_commitment.get_lde_values(i, step);
                let local_zs_partial_and_lookup =
                    zs_partial_products_and_lookup_commitment.get_lde_values(i, step);
                let next_zs_partial_and_lookup =
                    zs_partial_products_and_lookup_commitment.get_lde_values(i_next, step);

                let local_zs = &local_zs_partial_and_lookup[common_data.zs_range()];

                let next_zs = &next_zs_partial_and_lookup[common_data.zs_range()];

                let partial_products =
                    &local_zs_partial_and_lookup[common_data.partial_products_range()];

                if has_lookup {
                    let local_lookup_zs = &local_zs_partial_and_lookup[common_data.lookup_range()];

                    let next_lookup_zs = &next_zs_partial_and_lookup[common_data.lookup_range()];
                    debug_assert_eq!(local_lookup_zs.len(), common_data.num_all_lookup_polys());

                    local_lookup_batch.push(local_lookup_zs);
                    next_lookup_batch.push(next_lookup_zs);
                }

                debug_assert_eq!(local_wires.len(), common_data.config.num_wires);
                debug_assert_eq!(local_zs.len(), num_challenges);

                local_constants_batch_refs.push(local_constants);
                local_wires_batch_refs.push(local_wires);

                shifted_xs_batch.push(shifted_x);
                local_zs_batch.push(local_zs);
                next_zs_batch.push(next_zs);
                partial_products_batch.push(partial_products);
                s_sigmas_batch.push(s_sigmas);
            }

            // NB (JN): I'm not sure how (in)efficient the below is. It needs measuring.
            let mut local_constants_batch =
                vec![F::ZERO; xs_batch.len() * local_constants_batch_refs[0].len()];
            for i in 0..local_constants_batch_refs[0].len() {
                for (j, constants) in local_constants_batch_refs.iter().enumerate() {
                    local_constants_batch[i * xs_batch.len() + j] = constants[i];
                }
            }

            let mut local_wires_batch =
                vec![F::ZERO; xs_batch.len() * local_wires_batch_refs[0].len()];
            for i in 0..local_wires_batch_refs[0].len() {
                for (j, wires) in local_wires_batch_refs.iter().enumerate() {
                    local_wires_batch[i * xs_batch.len() + j] = wires[i];
                }
            }

            let vars_batch = EvaluationVarsBaseBatch::new(
                xs_batch.len(),
                &local_constants_batch,
                &local_wires_batch,
                public_inputs_hash,
            );

            let mut quotient_values_batch = eval_vanishing_poly_base_batch::<F, D>(
                common_data,
                &indices_batch,
                &shifted_xs_batch,
                vars_batch,
                &local_zs_batch,
                &next_zs_batch,
                &local_lookup_batch,
                &next_lookup_batch,
                &partial_products_batch,
                &s_sigmas_batch,
                betas,
                gammas,
                deltas,
                alphas,
                &z_h_on_coset,
            );

            for (&i, quotient_values) in indices_batch.iter().zip(quotient_values_batch.iter_mut())
            {
                let denominator_inv = z_h_on_coset.eval_inverse(i);
                quotient_values
                    .iter_mut()
                    .for_each(|v| *v *= denominator_inv);
            }
            quotient_values_batch
        })
        .collect();

    transpose(&quotient_values)
        .into_par_iter()
        .map(PolynomialValues::new)
        .map(|values| values.coset_ifft(F::coset_shift()))
        .collect()
}
```

`eval_vanishing_poly_base_batch` 也是个复杂的函数，它先调用 `evaluate_gate_constraints_base_batch` 计算各个 `gate` 带来的约束（`gate` 是否正确计算）在 batch 中的取值，再按 batch 计算其他约束：`vanishing_z1_terms`， `	vanishing_partial_products_terms`（调用 `check_partial_products` 计算），`vanishing_all_lookup_terms`（本例中没有）.

`vanishing_z1_terms` 对应的约束为 $Z(h)=1$，`vanishing_partial_products_terms` 对应的约束为 $\pi(x)$ 之间的递推关系.

最后调用 `reduce_with_powers_multi` 将它们组装起来，返回约束在多个点上的取值.

```rust
/// Like `eval_vanishing_poly`, but specialized for base field points. Batched.
pub(crate) fn eval_vanishing_poly_base_batch<F: RichField + Extendable<D>, const D: usize>(
    common_data: &CommonCircuitData<F, D>,
    indices_batch: &[usize],
    xs_batch: &[F],
    vars_batch: EvaluationVarsBaseBatch<F>,
    local_zs_batch: &[&[F]],
    next_zs_batch: &[&[F]],
    local_lookup_zs_batch: &[&[F]],
    next_lookup_zs_batch: &[&[F]],
    partial_products_batch: &[&[F]],
    s_sigmas_batch: &[&[F]],
    betas: &[F],
    gammas: &[F],
    deltas: &[F],
    alphas: &[F],
    z_h_on_coset: &ZeroPolyOnCoset<F>,
) -> Vec<Vec<F>> {
    let has_lookup = common_data.num_lookup_polys != 0;

    let n = indices_batch.len();
    assert_eq!(xs_batch.len(), n);
    assert_eq!(vars_batch.len(), n);
    assert_eq!(local_zs_batch.len(), n);
    assert_eq!(next_zs_batch.len(), n);
    if has_lookup {
        assert_eq!(local_lookup_zs_batch.len(), n);
        assert_eq!(next_lookup_zs_batch.len(), n);
    } else {
        assert_eq!(local_lookup_zs_batch.len(), 0);
        assert_eq!(next_lookup_zs_batch.len(), 0);
    }
    assert_eq!(partial_products_batch.len(), n);
    assert_eq!(s_sigmas_batch.len(), n);

    let max_degree = common_data.quotient_degree_factor;
    let num_prods = common_data.num_partial_products;

    let num_gate_constraints = common_data.num_gate_constraints;

    let constraint_terms_batch =
        evaluate_gate_constraints_base_batch::<F, D>(common_data, vars_batch);
    debug_assert!(constraint_terms_batch.len() == n * num_gate_constraints);

    let num_challenges = common_data.config.num_challenges;
    let num_routed_wires = common_data.config.num_routed_wires;

    let mut numerator_values = Vec::with_capacity(num_routed_wires);
    let mut denominator_values = Vec::with_capacity(num_routed_wires);

    // The L_0(x) (Z(x) - 1) vanishing terms.
    let mut vanishing_z_1_terms = Vec::with_capacity(num_challenges);
    // The terms checking the partial products.
    let mut vanishing_partial_products_terms = Vec::new();

    // The terms checking the lookup constraints.
    let mut vanishing_all_lookup_terms = if has_lookup {
        let num_sldc_polys = common_data.num_lookup_polys - 1;
        Vec::with_capacity(
            common_data.config.num_challenges * (4 + common_data.luts.len() + 2 * num_sldc_polys),
        )
    } else {
        Vec::new()
    };

    let mut res_batch: Vec<Vec<F>> = Vec::with_capacity(n);
    for k in 0..n {
        let index = indices_batch[k];
        let x = xs_batch[k];
        let vars = vars_batch.view(k);

        let lookup_selectors: Vec<F> = (0..common_data.num_lookup_selectors)
            .map(|i| vars.local_constants[common_data.selectors_info.num_selectors() + i])
            .collect();

        let local_zs = local_zs_batch[k];
        let next_zs = next_zs_batch[k];
        let local_lookup_zs = if has_lookup {
            local_lookup_zs_batch[k]
        } else {
            &[]
        };

        let next_lookup_zs = if has_lookup {
            next_lookup_zs_batch[k]
        } else {
            &[]
        };

        let partial_products = partial_products_batch[k];
        let s_sigmas = s_sigmas_batch[k];

        let constraint_terms = PackedStridedView::new(&constraint_terms_batch, n, k);

        let l_0_x = z_h_on_coset.eval_l_0(index, x);
        for i in 0..num_challenges {
            let z_x = local_zs[i];
            let z_gx = next_zs[i];
            vanishing_z_1_terms.push(l_0_x * z_x.sub_one());

            // If there are lookups in the circuit, then we add the lookup constraints.
            if has_lookup {
                let cur_deltas = &deltas[NUM_COINS_LOOKUP * i..NUM_COINS_LOOKUP * (i + 1)];

                let cur_local_lookup_zs = &local_lookup_zs
                    [common_data.num_lookup_polys * i..common_data.num_lookup_polys * (i + 1)];
                let cur_next_lookup_zs = &next_lookup_zs
                    [common_data.num_lookup_polys * i..common_data.num_lookup_polys * (i + 1)];

                let lookup_constraints = check_lookup_constraints_batch(
                    common_data,
                    vars,
                    cur_local_lookup_zs,
                    cur_next_lookup_zs,
                    &lookup_selectors,
                    cur_deltas.try_into().unwrap(),
                );
                vanishing_all_lookup_terms.extend(lookup_constraints);
            }

            numerator_values.extend((0..num_routed_wires).map(|j| {
                let wire_value = vars.local_wires[j];
                let k_i = common_data.k_is[j];
                let s_id = k_i * x;
                wire_value + betas[i] * s_id + gammas[i]
            }));
            denominator_values.extend((0..num_routed_wires).map(|j| {
                let wire_value = vars.local_wires[j];
                let s_sigma = s_sigmas[j];
                wire_value + betas[i] * s_sigma + gammas[i]
            }));

            // The partial products considered for this iteration of `i`.
            let current_partial_products = &partial_products[i * num_prods..(i + 1) * num_prods];
            // Check the numerator partial products.
            let partial_product_checks = check_partial_products(
                &numerator_values,
                &denominator_values,
                current_partial_products,
                z_x,
                z_gx,
                max_degree,
            );
            vanishing_partial_products_terms.extend(partial_product_checks);

            numerator_values.clear();
            denominator_values.clear();
        }

        let vanishing_terms = vanishing_z_1_terms
            .iter()
            .chain(vanishing_partial_products_terms.iter())
            .chain(vanishing_all_lookup_terms.iter())
            .chain(constraint_terms);
        let res = plonk_common::reduce_with_powers_multi(vanishing_terms, alphas);
        res_batch.push(res);

        vanishing_z_1_terms.clear();
        vanishing_partial_products_terms.clear();
        vanishing_all_lookup_terms.clear();
    }
    res_batch
}
```

`evaluate_gate_constraints_base_batch` 对每个门调用 `eval_filtered_base_batch`，计算得到当前 batch 相应的约束，并将其加入 `gate_constraints_batch` 中（加入时调用 `batch_add_inplace` 叠加，最后总约束个数为与最多的门相同）

```rust
/// Evaluate all gate constraints in the base field.
///
/// Returns a vector of `num_gate_constraints * vars_batch.len()` field elements. The constraints
/// corresponding to `vars_batch[i]` are found in `result[i], result[vars_batch.len() + i],
/// result[2 * vars_batch.len() + i], ...`.
pub fn evaluate_gate_constraints_base_batch<F: RichField + Extendable<D>, const D: usize>(
    common_data: &CommonCircuitData<F, D>,
    vars_batch: EvaluationVarsBaseBatch<F>,
) -> Vec<F> {
    let mut constraints_batch = vec![F::ZERO; common_data.num_gate_constraints * vars_batch.len()];
    for (i, gate) in common_data.gates.iter().enumerate() {
        let selector_index = common_data.selectors_info.selector_indices[i];
        let gate_constraints_batch = gate.0.eval_filtered_base_batch(
            vars_batch,
            i,
            selector_index,
            common_data.selectors_info.groups[selector_index].clone(),
            common_data.selectors_info.num_selectors(),
            common_data.num_lookup_selectors,
        );
        debug_assert!(
            gate_constraints_batch.len() <= constraints_batch.len(),
            "num_constraints() gave too low of a number"
        );
        // below adds all constraints for all points
        batch_add_inplace(
            &mut constraints_batch[..gate_constraints_batch.len()],
            &gate_constraints_batch,
        );
    }
    constraints_batch
}
```

`eval_filtered_base_batch` 会先调用 `compute_filter` 计算 batch 中每个点对应的 filter，再乘以当前门的约束. 若该点对应当前门，则约束算式应为 0，若该点不对应当前门，则 filter 应为 0.

实际调试时发现绝大多数约束的计算结果都不为 0，原因==有待研究==.

```rust
/// The result is an array of length `vars_batch.len() * self.num_constraints()`. Constraint `j`
/// for point `i` is at index `j * batch_size + i`.
fn eval_filtered_base_batch(
    &self,
    mut vars_batch: EvaluationVarsBaseBatch<F>,
    row: usize,
    selector_index: usize,
    group_range: Range<usize>,
    num_selectors: usize,
    num_lookup_selectors: usize,
) -> Vec<F> {
    let filters: Vec<_> = vars_batch
        .iter()
        .map(|vars| {
            compute_filter(
                row,
                group_range.clone(),
                vars.local_constants[selector_index],
                num_selectors > 1,
            )
        })
        .collect();
    vars_batch.remove_prefix(num_selectors + num_lookup_selectors);
    let mut res_batch = self.eval_unfiltered_base_batch(vars_batch);
    for res_chunk in res_batch.chunks_exact_mut(filters.len()) {
        batch_multiply_inplace(res_chunk, &filters);
    }
    res_batch
}
```

`compute_filter` 的计算为 $({\rm{UNUSED}}-s)\underset{i=0\\i\neq{k}}{\overset{n-1}{\Pi}}(i-s)$，其中 k 为当前行，s 为输入的点:

```rust
/// A gate's filter designed so that it is non-zero if `s = row`.
fn compute_filter<K: Field>(row: usize, group_range: Range<usize>, s: K, many_selector: bool) -> K {
    debug_assert!(group_range.contains(&row));
    group_range
        .filter(|&i| i != row)
        .chain(many_selector.then_some(UNUSED_SELECTOR))
        .map(|i| K::from_canonical_usize(i) - s)
        .product()
}
```

`eval_unfiltered_base_batch` 会调用 `self.eval_unfiltered_base_packed` 计算当前门的约束，返回一个数组，每 `vars_batch.len()` 个数对应同一个约束:

```rust
fn eval_unfiltered_base_batch(&self, vars_base: EvaluationVarsBaseBatch<F>) -> Vec<F> {
    self.eval_unfiltered_base_batch_packed(vars_base)
}

/// Evaluates entire batch of points. Returns a matrix of constraints. Constraint `j` for point
/// `i` is at `index j * batch_size + i`.
fn eval_unfiltered_base_batch_packed(&self, vars_batch: EvaluationVarsBaseBatch<F>) -> Vec<F> {
    let mut res = vec![F::ZERO; vars_batch.len() * self.num_constraints()];
    let (vars_packed_iter, vars_leftovers_iter) = vars_batch.pack::<<F as Packable>::Packing>();
    let leftovers_start = vars_batch.len() - vars_leftovers_iter.len();
    for (i, vars_packed) in vars_packed_iter.enumerate() {
        self.eval_unfiltered_base_packed(
            vars_packed,
            StridedConstraintConsumer::new(
                &mut res[..],
                vars_batch.len(),
                <F as Packable>::Packing::WIDTH * i,
            ),
        );
    }
    for (i, vars_leftovers) in vars_leftovers_iter.enumerate() {
        self.eval_unfiltered_base_packed(
            vars_leftovers,
            StridedConstraintConsumer::new(&mut res[..], vars_batch.len(), leftovers_start + i),
        );
    }
    res
}

fn eval_unfiltered_base_packed<P: PackedField<Scalar = F>>(
    &self,
    vars_base: EvaluationVarsBasePacked<P>,
    yield_constr: StridedConstraintConsumer<P>,
);
```

`eval_unfiltered_base_packed` 对不同门有不同实现，都是在对 commitment 中的值运算后相减，表示它们满足相应约束:

```rust
impl<F: RichField + Extendable<D>, const D: usize> PackedEvaluableBase<F, D> for ConstantGate {
    fn eval_unfiltered_base_packed<P: PackedField<Scalar = F>>(
        &self,
        vars: EvaluationVarsBasePacked<P>,
        mut yield_constr: StridedConstraintConsumer<P>,
    ) {
        yield_constr.many((0..self.num_consts).map(|i| {
            vars.local_constants[self.const_input(i)] - vars.local_wires[self.wire_output(i)]
        }));
    }
}

impl<F: RichField + Extendable<D>, const D: usize> PackedEvaluableBase<F, D> for PublicInputGate {
    fn eval_unfiltered_base_packed<P: PackedField<Scalar = F>>(
        &self,
        vars: EvaluationVarsBasePacked<P>,
        mut yield_constr: StridedConstraintConsumer<P>,
    ) {
        yield_constr.many(
            Self::wires_public_inputs_hash()
                .zip(vars.public_inputs_hash.elements)
                .map(|(wire, hash_part)| vars.local_wires[wire] - hash_part),
        )
    }
}

impl<F: RichField + Extendable<D>, const D: usize> PackedEvaluableBase<F, D> for ArithmeticGate {
    fn eval_unfiltered_base_packed<P: PackedField<Scalar = F>>(
        &self,
        vars: EvaluationVarsBasePacked<P>,
        mut yield_constr: StridedConstraintConsumer<P>,
    ) {
        let const_0 = vars.local_constants[0];
        let const_1 = vars.local_constants[1];

        for i in 0..self.num_ops {
            let multiplicand_0 = vars.local_wires[Self::wire_ith_multiplicand_0(i)];
            let multiplicand_1 = vars.local_wires[Self::wire_ith_multiplicand_1(i)];
            let addend = vars.local_wires[Self::wire_ith_addend(i)];
            let output = vars.local_wires[Self::wire_ith_output(i)];
            let computed_output = multiplicand_0 * multiplicand_1 * const_0 + addend * const_1;

            yield_constr.one(output - computed_output);
        }
    }
}
```

`check_partial_products` 通过减法来计算 $Z(x),\pi_i(x),Z(gx)$ 之间的递推关系:

```rust
/// Checks the relationship between each pair of partial product accumulators. In particular, this
/// sequence of accumulators starts with `Z(x)`, then contains each partial product polynomials
/// `p_i(x)`, and finally `Z(g x)`. See the partial products section of the Plonky2 paper.
pub(crate) fn check_partial_products<F: Field>(
    numerators: &[F],
    denominators: &[F],
    partials: &[F],
    z_x: F,
    z_gx: F,
    max_degree: usize,
) -> Vec<F> {
    debug_assert!(max_degree > 1);
    let product_accs = iter::once(&z_x)
        .chain(partials.iter())
        .chain(iter::once(&z_gx));
    let chunk_size = max_degree;
    numerators
        .chunks(chunk_size)
        .zip_eq(denominators.chunks(chunk_size))
        .zip_eq(product_accs.tuple_windows())
        .map(|((nume_chunk, deno_chunk), (&prev_acc, &next_acc))| {
            let num_chunk_product = nume_chunk.iter().copied().product();
            let den_chunk_product = deno_chunk.iter().copied().product();
            // Assert that next_acc * deno_product = prev_acc * nume_product.
            prev_acc * num_chunk_product - next_acc * den_chunk_product
        })
        .collect()
}
```

#### Openings and Proof

生成 `OpeningSet` 以及相应的 `OpeningProof`.

首先从 `challenger` 获取打开的位置 `zeta` 以及 `g * zeta`（这里的 `g` 是单位根而不是生成元），再调用 `OpeningSet::new` 生成 `openings`，从 `challenger` 获取 FRI 的 `instance`，再调用 `PolynomialBatch::prove_openings` 生成 FRI 的 `opening_proof` :

```rust
challenger.observe_cap::<C::Hasher>(&quotient_polys_commitment.merkle_tree.cap);

let zeta = challenger.get_extension_challenge::<D>();
// To avoid leaking witness data, we want to ensure that our opening locations, `zeta` and
// `g * zeta`, are not in our subgroup `H`. It suffices to check `zeta` only, since
// `(g * zeta)^n = zeta^n`, where `n` is the order of `g`.
let g = F::Extension::primitive_root_of_unity(common_data.degree_bits());
ensure!(
    zeta.exp_power_of_2(common_data.degree_bits()) != F::Extension::ONE,
    "Opening point is in the subgroup."
);

let openings = timed!(
    timing,
    "construct the opening set, including lookups",
    OpeningSet::new(
        zeta,
        g,
        &prover_data.constants_sigmas_commitment,
        &wires_commitment,
        &partial_products_zs_and_lookup_commitment,
        &quotient_polys_commitment,
        common_data
    )
);
challenger.observe_openings(&openings.to_fri_openings());
let instance = common_data.get_fri_instance(zeta);

let opening_proof = timed!(
    timing,
    "compute opening proofs",
    PolynomialBatch::<F, C, D>::prove_openings(
        &instance,
        &[
            &prover_data.constants_sigmas_commitment,
            &wires_commitment,
            &partial_products_zs_and_lookup_commitment,
            &quotient_polys_commitment,
        ],
        &mut challenger,
        &common_data.fri_params,
        timing,
    )
);
```

最后，将几部分组装为最后的 `proof`，结束 `prove` 流程:

```rust
let proof = Proof::<F, C, D> {
    wires_cap: wires_commitment.merkle_tree.cap,
    plonk_zs_partial_products_cap: partial_products_zs_and_lookup_commitment.merkle_tree.cap,
    quotient_polys_cap: quotient_polys_commitment.merkle_tree.cap,
    openings,
    opening_proof,
};
Ok(ProofWithPublicInputs::<F, C, D> {
    proof,
    public_inputs,
})
```

##### Openings

在 `OpeningSet::new` 中，先定义了函数 `eval_commitment`，它会调用 `eval` 求出一个 `PolynomialBatch` 中每个多项式在给定点的值.

```rust
pub fn new<C: GenericConfig<D, F = F>>(
    zeta: F::Extension,
    g: F::Extension,
    constants_sigmas_commitment: &PolynomialBatch<F, C, D>,
    wires_commitment: &PolynomialBatch<F, C, D>,
    zs_partial_products_lookup_commitment: &PolynomialBatch<F, C, D>,
    quotient_polys_commitment: &PolynomialBatch<F, C, D>,
    common_data: &CommonCircuitData<F, D>,
) -> Self {
    let eval_commitment = |z: F::Extension, c: &PolynomialBatch<F, C, D>| {
        c.polynomials
        .par_iter()
        .map(|p| p.to_extension().eval(z))
        .collect::<Vec<_>>()
    };
    let constants_sigmas_eval = eval_commitment(zeta, constants_sigmas_commitment);

    // `zs_partial_products_lookup_eval` contains the permutation argument polynomials as well as lookup polynomials.
    let zs_partial_products_lookup_eval =
    eval_commitment(zeta, zs_partial_products_lookup_commitment);
    let zs_partial_products_lookup_next_eval =
    eval_commitment(g * zeta, zs_partial_products_lookup_commitment);
    let quotient_polys = eval_commitment(zeta, quotient_polys_commitment);

    Self {
        constants: constants_sigmas_eval[common_data.constants_range()].to_vec(),
        plonk_sigmas: constants_sigmas_eval[common_data.sigmas_range()].to_vec(),
        wires: eval_commitment(zeta, wires_commitment),
        plonk_zs: zs_partial_products_lookup_eval[common_data.zs_range()].to_vec(),
        plonk_zs_next: zs_partial_products_lookup_next_eval[common_data.zs_range()].to_vec(),
        partial_products: zs_partial_products_lookup_eval[common_data.partial_products_range()]
        	.to_vec(),
        quotient_polys,
        lookup_zs: zs_partial_products_lookup_eval[common_data.lookup_range()].to_vec(),
        lookup_zs_next: zs_partial_products_lookup_next_eval[common_data.lookup_range()]
        	.to_vec(),
    }
}
```

##### Instance

调用 `commondata.get_fri_instance` 得到 `instance` :

调用 `self.fri_all_polys`，`self.fri_next_batch_polys`，`self.fri_oracles` 得到一系列多项式信息和 `oracles` 信息:

```rust
pub(crate) fn get_fri_instance(&self, zeta: F::Extension) -> FriInstanceInfo<F, D> {
    // All polynomials are opened at zeta.
    let zeta_batch = FriBatchInfo {
        point: zeta,
        polynomials: self.fri_all_polys(),
    };

    // The Z polynomials are also opened at g * zeta.
    let g = F::Extension::primitive_root_of_unity(self.degree_bits());
    let zeta_next = g * zeta;
    let zeta_next_batch = FriBatchInfo {
        point: zeta_next,
        polynomials: self.fri_next_batch_polys(),
    };

    let openings = vec![zeta_batch, zeta_next_batch];
    FriInstanceInfo {
        oracles: self.fri_oracles(),
        batches: openings,
    }
}
```

`fri_all_polys` 一类函数都只是求出多项式或 `oracle` 的信息，并不涉及具体数值:

```rust
fn fri_all_polys(&self) -> Vec<FriPolynomialInfo> {
    [
        self.fri_preprocessed_polys(),
        self.fri_wire_polys(),
        self.fri_zs_partial_products_polys(),
        self.fri_quotient_polys(),
        self.fri_lookup_polys(),
    ]
    .concat()
}

fn fri_preprocessed_polys(&self) -> Vec<FriPolynomialInfo> {
    FriPolynomialInfo::from_range(
        PlonkOracle::CONSTANTS_SIGMAS.index,
        0..self.num_preprocessed_polys(),
    )
}

/// Returns polynomials that require evaluation at `zeta` and `g * zeta`.
fn fri_next_batch_polys(&self) -> Vec<FriPolynomialInfo> {
    [self.fri_zs_polys(), self.fri_lookup_polys()].concat()
}

fn fri_oracles(&self) -> Vec<FriOracleInfo> {
    vec![
        FriOracleInfo {
            num_polys: self.num_preprocessed_polys(),
            blinding: PlonkOracle::CONSTANTS_SIGMAS.blinding,
        },
        FriOracleInfo {
            num_polys: self.config.num_wires,
            blinding: PlonkOracle::WIRES.blinding,
        },
        FriOracleInfo {
            num_polys: self.num_zs_partial_products_polys() + self.num_all_lookup_polys(),
            blinding: PlonkOracle::ZS_PARTIAL_PRODUCTS.blinding,
        },
        FriOracleInfo {
            num_polys: self.num_quotient_polys(),
            blinding: PlonkOracle::QUOTIENT.blinding,
        },
    ]
}
```

##### FRI Proof

`prove_openings` 生成随机数 `alpha`，再类似之前方法调用 `alpha.reduce_polys_base` 计算 FRI 的 `final_poly` 系数，求出它的低度扩展，最后调用 `fri_proof` 生成 `proof` :

```rust
/// Produces a batch opening proof.
pub fn prove_openings(
    instance: &FriInstanceInfo<F, D>,
    oracles: &[&Self],
    challenger: &mut Challenger<F, C::Hasher>,
    fri_params: &FriParams,
    timing: &mut TimingTree,
) -> FriProof<F, C::Hasher, D> {
    assert!(D > 1, "Not implemented for D=1.");
    let alpha = challenger.get_extension_challenge::<D>();
    let mut alpha = ReducingFactor::new(alpha);

    // Final low-degree polynomial that goes into FRI.
    let mut final_poly = PolynomialCoeffs::empty();

    // Each batch `i` consists of an opening point `z_i` and polynomials `{f_ij}_j` to be opened at that point.
    // For each batch, we compute the composition polynomial `F_i = sum alpha^j f_ij`,
    // where `alpha` is a random challenge in the extension field.
    // The final polynomial is then computed as `final_poly = sum_i alpha^(k_i) (F_i(X) - F_i(z_i))/(X-z_i)`
    // where the `k_i`s are chosen such that each power of `alpha` appears only once in the final sum.
    // There are usually two batches for the openings at `zeta` and `g * zeta`.
    // The oracles used in Plonky2 are given in `FRI_ORACLES` in `plonky2/src/plonk/plonk_common.rs`.
    for FriBatchInfo { point, polynomials } in &instance.batches {
        // Collect the coefficients of all the polynomials in `polynomials`.
        let polys_coeff = polynomials.iter().map(|fri_poly| {
            &oracles[fri_poly.oracle_index].polynomials[fri_poly.polynomial_index]
        });
        let composition_poly = timed!(
            timing,
            &format!("reduce batch of {} polynomials", polynomials.len()),
            alpha.reduce_polys_base(polys_coeff)
        );
        let mut quotient = composition_poly.divide_by_linear(*point);
        quotient.coeffs.push(F::Extension::ZERO); // pad back to power of two
        alpha.shift_poly(&mut final_poly);
        final_poly += quotient;
    }

    let lde_final_poly = final_poly.lde(fri_params.config.rate_bits);
    let lde_final_values = timed!(
        timing,
        &format!("perform final FFT {}", lde_final_poly.len()),
        lde_final_poly.coset_fft(F::coset_shift().into())
    );

    let fri_proof = fri_proof::<F, C, D>(
        &oracles
        .par_iter()
        .map(|c| &c.merkle_tree)
        .collect::<Vec<_>>(),
        lde_final_poly,
        lde_final_values,
        challenger,
        fri_params,
        timing,
    );

    fri_proof
}
```

`reduce_polys_base` 会将一个 batch 的多项式系数乘上 `alpha` 的幂次，再加起来得到 `composition_poly`，即 $F_i(x)=\underset{j}\sum\alpha^jf_{ij}(x)$:

之后再调用 `divide_by_poly` 计算得到 $(F_i(x)-F_i(z_i))/(x-z_i)$，组装成 `final_poly` 即 $\underset{i}\sum\alpha^{k_i}(F_i(x)-F_i(z_i))/(x-z_i)$.

```rust
pub fn reduce_polys_base<BF: Extendable<D, Extension = F>, const D: usize>(
    &mut self,
    polys: impl IntoIterator<Item = impl Borrow<PolynomialCoeffs<BF>>>,
) -> PolynomialCoeffs<F> {
    self.base
        .powers()
        .zip(polys)
        .map(|(base_power, poly)| {
            self.count += 1;
            poly.borrow().mul_extension(base_power)
        })
        .sum()
}
```

`fri_proof` 包括三个阶段: Commit phase, Proof of Work phase, Query phase.

分别调用 `fri_committed_trees`, `fri_proof_of_work`, `fri_prover_query_rounds` 并将结果组装为 `FriProof`.

```rust
/// Builds a FRI proof.
pub fn fri_proof<F: RichField + Extendable<D>, C: GenericConfig<D, F = F>, const D: usize>(
    initial_merkle_trees: &[&MerkleTree<F, C::Hasher>],
    // Coefficients of the polynomial on which the LDT is performed. Only the first `1/rate` coefficients are non-zero.
    lde_polynomial_coeffs: PolynomialCoeffs<F::Extension>,
    // Evaluation of the polynomial on the large domain.
    lde_polynomial_values: PolynomialValues<F::Extension>,
    challenger: &mut Challenger<F, C::Hasher>,
    fri_params: &FriParams,
    timing: &mut TimingTree,
) -> FriProof<F, C::Hasher, D> {
    let n = lde_polynomial_values.len();
    assert_eq!(lde_polynomial_coeffs.len(), n);

    // Commit phase
    let (trees, final_coeffs) = timed!(
        timing,
        "fold codewords in the commitment phase",
        fri_committed_trees::<F, C, D>(
            lde_polynomial_coeffs,
            lde_polynomial_values,
            challenger,
            fri_params,
        )
    );

    // PoW phase
    let pow_witness = timed!(
        timing,
        "find proof-of-work witness",
        fri_proof_of_work::<F, C, D>(challenger, &fri_params.config)
    );

    // Query phase
    let query_round_proofs =
        fri_prover_query_rounds::<F, C, D>(initial_merkle_trees, &trees, challenger, n, fri_params);

    FriProof {
        commit_phase_merkle_caps: trees.iter().map(|t| t.cap.clone()).collect(),
        query_round_proofs,
        final_poly: final_coeffs,
        pow_witness,
    }
}
```

在 `fri_committed_trees` 中，对 `fri_params.reduction_arity_bits` 中的每一个 `arity_bits`（一般为 1），调用 `flatten`

```rust
fn fri_committed_trees<F: RichField + Extendable<D>, C: GenericConfig<D, F = F>, const D: usize>(
    mut coeffs: PolynomialCoeffs<F::Extension>,
    mut values: PolynomialValues<F::Extension>,
    challenger: &mut Challenger<F, C::Hasher>,
    fri_params: &FriParams,
) -> FriCommitedTrees<F, C, D> {
    let mut trees = Vec::new();

    let mut shift = F::MULTIPLICATIVE_GROUP_GENERATOR;
    for arity_bits in &fri_params.reduction_arity_bits {
        let arity = 1 << arity_bits;

        reverse_index_bits_in_place(&mut values.values);
        let chunked_values = values
            .values
            .par_chunks(arity)
            .map(|chunk: &[F::Extension]| flatten(chunk))
            .collect();
        let tree = MerkleTree::<F, C::Hasher>::new(chunked_values, fri_params.config.cap_height);

        challenger.observe_cap(&tree.cap);
        trees.push(tree);

        let beta = challenger.get_extension_challenge::<D>();
        // P(x) = sum_{i<r} x^i * P_i(x^r) becomes sum_{i<r} beta^i * P_i(x).
        coeffs = PolynomialCoeffs::new(
            coeffs
                .coeffs
                .par_chunks_exact(arity)
                .map(|chunk| reduce_with_powers(chunk, beta))
                .collect::<Vec<_>>(),
        );
        shift = shift.exp_u64(arity as u64);
        values = coeffs.coset_fft(shift.into())
    }

    // The coefficients being removed here should always be zero.
    coeffs
        .coeffs
        .truncate(coeffs.len() >> fri_params.config.rate_bits);

    challenger.observe_extension_elements(&coeffs.coeffs);
    (trees, coeffs)
}
```



```rust
/// Performs the proof-of-work (a.k.a. grinding) step of the FRI protocol. Returns the PoW witness.
fn fri_proof_of_work<F: RichField + Extendable<D>, C: GenericConfig<D, F = F>, const D: usize>(
    challenger: &mut Challenger<F, C::Hasher>,
    config: &FriConfig,
) -> F {
    let min_leading_zeros = config.proof_of_work_bits + (64 - F::order().bits()) as u32;

    // The easiest implementation would be repeatedly clone our Challenger. With each clone, we'd
    // observe an incrementing PoW witness, then get the PoW response. If it contained sufficient
    // leading zeros, we'd end the search, and store this clone as our new challenger.
    //
    // However, performance is critical here. We want to avoid cloning Challenger, particularly
    // since it stores vectors, which means allocations. We'd like a more compact state to clone.
    //
    // We know that a duplex will be performed right after we send the PoW witness, so we can ignore
    // any output_buffer, which will be invalidated. We also know
    // input_buffer.len() < H::Permutation::WIDTH, an invariant of Challenger.
    //
    // We separate the duplex operation into two steps, one which can be performed now, and the
    // other which depends on the PoW witness candidate. The first step is the overwrite our sponge
    // state with any inputs (excluding the PoW witness candidate). The second step is to overwrite
    // one more element of our sponge state with the candidate, then apply the permutation,
    // obtaining our duplex's post-state which contains the PoW response.
    let mut duplex_intermediate_state = challenger.sponge_state;
    let witness_input_pos = challenger.input_buffer.len();
    duplex_intermediate_state.set_from_iter(challenger.input_buffer.clone().into_iter(), 0);

    let pow_witness = (0..=F::NEG_ONE.to_canonical_u64())
        .into_par_iter()
        .find_any(|&candidate| {
            let mut duplex_state = duplex_intermediate_state;
            duplex_state.set_elt(F::from_canonical_u64(candidate), witness_input_pos);
            duplex_state.permute();
            let pow_response = duplex_state.squeeze().iter().last().unwrap();
            let leading_zeros = pow_response.to_canonical_u64().leading_zeros();
            leading_zeros >= min_leading_zeros
        })
        .map(F::from_canonical_u64)
        .expect("Proof of work failed. This is highly unlikely!");

    // Recompute pow_response using our normal Challenger code, and make sure it matches.
    challenger.observe_element(pow_witness);
    let pow_response = challenger.get_challenge();
    let leading_zeros = pow_response.to_canonical_u64().leading_zeros();
    assert!(leading_zeros >= min_leading_zeros);
    pow_witness
}
```



```rust
fn fri_prover_query_rounds<
    F: RichField + Extendable<D>,
    C: GenericConfig<D, F = F>,
    const D: usize,
>(
    initial_merkle_trees: &[&MerkleTree<F, C::Hasher>],
    trees: &[MerkleTree<F, C::Hasher>],
    challenger: &mut Challenger<F, C::Hasher>,
    n: usize,
    fri_params: &FriParams,
) -> Vec<FriQueryRound<F, C::Hasher, D>> {
    challenger
        .get_n_challenges(fri_params.config.num_query_rounds)
        .into_par_iter()
        .map(|rand| {
            let x_index = rand.to_canonical_u64() as usize % n;
            fri_prover_query_round::<F, C, D>(initial_merkle_trees, trees, x_index, fri_params)
        })
        .collect()
}

fn fri_prover_query_round<
    F: RichField + Extendable<D>,
    C: GenericConfig<D, F = F>,
    const D: usize,
>(
    initial_merkle_trees: &[&MerkleTree<F, C::Hasher>],
    trees: &[MerkleTree<F, C::Hasher>],
    mut x_index: usize,
    fri_params: &FriParams,
) -> FriQueryRound<F, C::Hasher, D> {
    let mut query_steps = Vec::new();
    let initial_proof = initial_merkle_trees
        .iter()
        .map(|t| (t.get(x_index).to_vec(), t.prove(x_index)))
        .collect::<Vec<_>>();
    for (i, tree) in trees.iter().enumerate() {
        let arity_bits = fri_params.reduction_arity_bits[i];
        let evals = unflatten(tree.get(x_index >> arity_bits));
        let merkle_proof = tree.prove(x_index >> arity_bits);

        query_steps.push(FriQueryStep {
            evals,
            merkle_proof,
        });

        x_index >>= arity_bits;
    }
    FriQueryRound {
        initial_trees_proof: FriInitialTreeProof {
            evals_proofs: initial_proof,
        },
        steps: query_steps,
    }
}
```









### Verify

通过调用 `data.verify(proof)` 进而调用 `verifier::verify`，我们利用 `VerifierOnlyCircuitData` 实现了对 `proof` 的验证:

`verify` 函数很简短，做了四件事：验证 `proof` 格式；计算哈希 `public_inputs_hash`；生成 `challenges`；用 `challenges` 验证.

其中 `verify_with_challenges` 最复杂.

```rust
pub(crate) fn verify<F: RichField + Extendable<D>, C: GenericConfig<D, F = F>, const D: usize>(
    proof_with_pis: ProofWithPublicInputs<F, C, D>,
    verifier_data: &VerifierOnlyCircuitData<C, D>,
    common_data: &CommonCircuitData<F, D>,
) -> Result<()> {
    validate_proof_with_pis_shape(&proof_with_pis, common_data)?;

    let public_inputs_hash = proof_with_pis.get_public_inputs_hash();
    let challenges = proof_with_pis.get_challenges(
        public_inputs_hash,
        &verifier_data.circuit_digest,
        common_data,
    )?;

    verify_with_challenges::<F, C, D>(
        proof_with_pis.proof,
        public_inputs_hash,
        challenges,
        verifier_data,
        common_data,
    )
}
```

#### Shape

`validate_proof_with_pis_shape` 调用 `validate_proof_shape` 验证各 merkle tree cap 的高度、各向量的长度等，再验证 `public_inputs` 的大小:

```rust
pub(crate) fn validate_proof_with_pis_shape<F, C, const D: usize>(
    proof_with_pis: &ProofWithPublicInputs<F, C, D>,
    common_data: &CommonCircuitData<F, D>,
) -> anyhow::Result<()>
where
    F: RichField + Extendable<D>,
    C: GenericConfig<D, F = F>,
{
    let ProofWithPublicInputs {
        proof,
        public_inputs,
    } = proof_with_pis;
    validate_proof_shape(proof, common_data)?;
    ensure!(
        public_inputs.len() == common_data.num_public_inputs,
        "Number of public inputs doesn't match circuit data."
    );
    Ok(())
}

fn validate_proof_shape<F, C, const D: usize>(
    proof: &Proof<F, C, D>,
    common_data: &CommonCircuitData<F, D>,
) -> anyhow::Result<()>
where
    F: RichField + Extendable<D>,
    C: GenericConfig<D, F = F>,
{
    let config = &common_data.config;
    let Proof {
        wires_cap,
        plonk_zs_partial_products_cap,
        quotient_polys_cap,
        openings,
        // The shape of the opening proof will be checked in the FRI verifier (see
        // validate_fri_proof_shape), so we ignore it here.
        opening_proof: _,
    } = proof;
    let OpeningSet {
        constants,
        plonk_sigmas,
        wires,
        plonk_zs,
        plonk_zs_next,
        partial_products,
        quotient_polys,
        lookup_zs,
        lookup_zs_next,
    } = openings;
    let cap_height = common_data.fri_params.config.cap_height;
    ensure!(wires_cap.height() == cap_height);
    ensure!(plonk_zs_partial_products_cap.height() == cap_height);
    ensure!(quotient_polys_cap.height() == cap_height);
    ensure!(constants.len() == common_data.num_constants);
    ensure!(plonk_sigmas.len() == config.num_routed_wires);
    ensure!(wires.len() == config.num_wires);
    ensure!(plonk_zs.len() == config.num_challenges);
    ensure!(plonk_zs_next.len() == config.num_challenges);
    ensure!(partial_products.len() == config.num_challenges * common_data.num_partial_products);
    ensure!(quotient_polys.len() == common_data.num_quotient_polys());
    ensure!(lookup_zs.len() == common_data.num_all_lookup_polys());
    ensure!(lookup_zs_next.len() == common_data.num_all_lookup_polys());
    Ok(())
}
```

#### Challenges

首先调用 `get_public_inputs_hash` 得到公共输入的 Poseidon hash:

```rust
pub fn get_public_inputs_hash(
    &self,
) -> <<C as GenericConfig<D>>::InnerHasher as Hasher<F>>::Hash {
    C::InnerHasher::hash_no_pad(&self.public_inputs)
}

fn hash_no_pad(input: &[F]) -> Self::Hash {
    hash_n_to_hash_no_pad::<F, Self::Permutation>(input)
}

pub fn hash_n_to_hash_no_pad<F: RichField, P: PlonkyPermutation<F>>(inputs: &[F]) -> HashOut<F> {
    HashOut::from_vec(hash_n_to_m_no_pad::<F, P>(inputs, NUM_HASH_OUT_ELTS))
}

/// Hash a message without any padding step. Note that this can enable length-extension attacks.
/// However, it is still collision-resistant in cases where the input has a fixed length.
pub fn hash_n_to_m_no_pad<F: RichField, P: PlonkyPermutation<F>>(
    inputs: &[F],
    num_outputs: usize,
) -> Vec<F> {
    let mut perm = P::new(repeat(F::ZERO));

    // Absorb all input chunks.
    for input_chunk in inputs.chunks(P::RATE) {
        perm.set_from_slice(input_chunk, 0);
        perm.permute();
    }

    // Squeeze until we have the desired number of outputs.
    let mut outputs = Vec::new();
    loop {
        for &item in perm.squeeze() {
            outputs.push(item);
            if outputs.len() == num_outputs {
                return outputs;
            }
        }
        perm.permute();
    }
}

fn permute(input: [Self; SPONGE_WIDTH]) -> [Self; SPONGE_WIDTH] {
    <F as Poseidon>::poseidon(input)
}
```

再调用 `get_challenges` 得到伪随机数，与 `prove` 中得到的一致:

| challenger.observe...         | challenger.get...            |
| ----------------------------- | ---------------------------- |
| circuit_digest                |                              |
| public_inputs_hash            |                              |
| wires_cap                     |                              |
|                               | plonk_betas                  |
|                               | plonk_gammas                 |
|                               | plonk_deltas (if has lookup) |
| plonk_zs_partial_products_cap |                              |
|                               | plonk_alphas                 |
| quotient_polys_cap            |                              |
|                               | plonk_zetas                  |
| openings.to_fri_openings()    |                              |
|                               | fri_alpha                    |
|                               | fri_betas                    |
| final_poly.coeffs             |                              |
| pow_witness                   |                              |
|                               | fri_pow_response             |
|                               | fri_query_indices            |

```rust
/// Computes all Fiat-Shamir challenges used in the Plonk proof.
pub fn get_challenges(
    &self,
    public_inputs_hash: <<C as GenericConfig<D>>::InnerHasher as Hasher<F>>::Hash,
    circuit_digest: &<<C as GenericConfig<D>>::Hasher as Hasher<C::F>>::Hash,
    common_data: &CommonCircuitData<F, D>,
) -> anyhow::Result<ProofChallenges<F, D>> {
    let Proof {
        wires_cap,
        plonk_zs_partial_products_cap,
        quotient_polys_cap,
        openings,
        opening_proof:
        FriProof {
            commit_phase_merkle_caps,
            final_poly,
            pow_witness,
            ..
        },
    } = &self.proof;

    get_challenges::<F, C, D>(
        public_inputs_hash,
        wires_cap,
        plonk_zs_partial_products_cap,
        quotient_polys_cap,
        openings,
        commit_phase_merkle_caps,
        final_poly,
        *pow_witness,
        circuit_digest,
        common_data,
    )
}

fn get_challenges<F: RichField + Extendable<D>, C: GenericConfig<D, F = F>, const D: usize>(
    public_inputs_hash: <<C as GenericConfig<D>>::InnerHasher as Hasher<F>>::Hash,
    wires_cap: &MerkleCap<F, C::Hasher>,
    plonk_zs_partial_products_cap: &MerkleCap<F, C::Hasher>,
    quotient_polys_cap: &MerkleCap<F, C::Hasher>,
    openings: &OpeningSet<F, D>,
    commit_phase_merkle_caps: &[MerkleCap<F, C::Hasher>],
    final_poly: &PolynomialCoeffs<F::Extension>,
    pow_witness: F,
    circuit_digest: &<<C as GenericConfig<D>>::Hasher as Hasher<C::F>>::Hash,
    common_data: &CommonCircuitData<F, D>,
) -> anyhow::Result<ProofChallenges<F, D>> {
    let config = &common_data.config;
    let num_challenges = config.num_challenges;

    let mut challenger = Challenger::<F, C::Hasher>::new();
    let has_lookup = common_data.num_lookup_polys != 0;

    // Observe the instance.
    challenger.observe_hash::<C::Hasher>(*circuit_digest);
    challenger.observe_hash::<C::InnerHasher>(public_inputs_hash);

    challenger.observe_cap::<C::Hasher>(wires_cap);
    let plonk_betas = challenger.get_n_challenges(num_challenges);
    let plonk_gammas = challenger.get_n_challenges(num_challenges);

    // If there are lookups in the circuit, we should get delta challenges as well.
    // But we can use the already generated `plonk_betas` and `plonk_gammas` as the first `plonk_deltas` challenges.
    let plonk_deltas = if has_lookup {
        let num_lookup_challenges = NUM_COINS_LOOKUP * num_challenges;
        let mut deltas = Vec::with_capacity(num_lookup_challenges);
        let num_additional_challenges = num_lookup_challenges - 2 * num_challenges;
        let additional = challenger.get_n_challenges(num_additional_challenges);
        deltas.extend(&plonk_betas);
        deltas.extend(&plonk_gammas);
        deltas.extend(additional);
        deltas
    } else {
        vec![]
    };

    // `plonk_zs_partial_products_cap` also contains the commitment to lookup polynomials.
    challenger.observe_cap::<C::Hasher>(plonk_zs_partial_products_cap);
    let plonk_alphas = challenger.get_n_challenges(num_challenges);

    challenger.observe_cap::<C::Hasher>(quotient_polys_cap);
    let plonk_zeta = challenger.get_extension_challenge::<D>();

    challenger.observe_openings(&openings.to_fri_openings());

    Ok(ProofChallenges {
        plonk_betas,
        plonk_gammas,
        plonk_alphas,
        plonk_deltas,
        plonk_zeta,
        fri_challenges: challenger.fri_challenges::<C, D>(
            commit_phase_merkle_caps,
            final_poly,
            pow_witness,
            common_data.degree_bits(),
            &config.fri_config,
        ),
    })
}

pub fn fri_challenges<C: GenericConfig<D, F = F>, const D: usize>(
    &mut self,
    commit_phase_merkle_caps: &[MerkleCap<F, C::Hasher>],
    final_poly: &PolynomialCoeffs<F::Extension>,
    pow_witness: F,
    degree_bits: usize,
    config: &FriConfig,
) -> FriChallenges<F, D>
	where
F: RichField + Extendable<D>,
{
    let num_fri_queries = config.num_query_rounds;
    let lde_size = 1 << (degree_bits + config.rate_bits);
    // Scaling factor to combine polynomials.
    let fri_alpha = self.get_extension_challenge::<D>();

    // Recover the random betas used in the FRI reductions.
    let fri_betas = commit_phase_merkle_caps
        .iter()
        .map(|cap| {
            self.observe_cap::<C::Hasher>(cap);
            self.get_extension_challenge::<D>()
        })
        .collect();

    self.observe_extension_elements(&final_poly.coeffs);

    self.observe_element(pow_witness);
    let fri_pow_response = self.get_challenge();

    let fri_query_indices = (0..num_fri_queries)
        .map(|_| self.get_challenge().to_canonical_u64() as usize % lde_size)
        .collect();

    FriChallenges {
        fri_alpha,
        fri_betas,
        fri_pow_response,
        fri_query_indices,
    }
}
```

#### Vanishing Polys

`verify_with_challenges` 调用 `eval_vanishing_poly` 计算 `vanishing_poly` 在 `zeta` 处的值，验证后再调用 `verify_fri_proof` :

验证时对每个 `vanishing_poly` 验证它在 $\zeta$ 处等于 $Z_H(\zeta)\cdot(t_0(\zeta)+t_1(\zeta)\cdot\zeta^n+t_2(\zeta)\cdot\zeta^{2n}+\cdots)$.

```rust
pub(crate) fn verify_with_challenges<
    F: RichField + Extendable<D>,
    C: GenericConfig<D, F = F>,
    const D: usize,
>(
    proof: Proof<F, C, D>,
    public_inputs_hash: <<C as GenericConfig<D>>::InnerHasher as Hasher<F>>::Hash,
    challenges: ProofChallenges<F, D>,
    verifier_data: &VerifierOnlyCircuitData<C, D>,
    common_data: &CommonCircuitData<F, D>,
) -> Result<()> {
    let local_constants = &proof.openings.constants;
    let local_wires = &proof.openings.wires;
    let vars = EvaluationVars {
        local_constants,
        local_wires,
        public_inputs_hash: &public_inputs_hash,
    };
    let local_zs = &proof.openings.plonk_zs;
    let next_zs = &proof.openings.plonk_zs_next;
    let local_lookup_zs = &proof.openings.lookup_zs;
    let next_lookup_zs = &proof.openings.lookup_zs_next;
    let s_sigmas = &proof.openings.plonk_sigmas;
    let partial_products = &proof.openings.partial_products;

    // Evaluate the vanishing polynomial at our challenge point, zeta.
    let vanishing_polys_zeta = eval_vanishing_poly::<F, D>(
        common_data,
        challenges.plonk_zeta,
        vars,
        local_zs,
        next_zs,
        local_lookup_zs,
        next_lookup_zs,
        partial_products,
        s_sigmas,
        &challenges.plonk_betas,
        &challenges.plonk_gammas,
        &challenges.plonk_alphas,
        &challenges.plonk_deltas,
    );

    // Check each polynomial identity, of the form `vanishing(x) = Z_H(x) quotient(x)`, at zeta.
    let quotient_polys_zeta = &proof.openings.quotient_polys;
    let zeta_pow_deg = challenges
        .plonk_zeta
        .exp_power_of_2(common_data.degree_bits());
    let z_h_zeta = zeta_pow_deg - F::Extension::ONE;
    // `quotient_polys_zeta` holds `num_challenges * quotient_degree_factor` evaluations.
    // Each chunk of `quotient_degree_factor` holds the evaluations of `t_0(zeta),...,t_{quotient_degree_factor-1}(zeta)`
    // where the "real" quotient polynomial is `t(X) = t_0(X) + t_1(X)*X^n + t_2(X)*X^{2n} + ...`.
    // So to reconstruct `t(zeta)` we can compute `reduce_with_powers(chunk, zeta^n)` for each
    // `quotient_degree_factor`-sized chunk of the original evaluations.
    for (i, chunk) in quotient_polys_zeta
        .chunks(common_data.quotient_degree_factor)
        .enumerate()
    {
        ensure!(vanishing_polys_zeta[i] == z_h_zeta * reduce_with_powers(chunk, zeta_pow_deg));
    }

    let merkle_caps = &[
        verifier_data.constants_sigmas_cap.clone(),
        proof.wires_cap,
        // In the lookup case, `plonk_zs_partial_products_cap` should also include the lookup commitment.
        proof.plonk_zs_partial_products_cap,
        proof.quotient_polys_cap,
    ];

    verify_fri_proof::<F, C, D>(
        &common_data.get_fri_instance(challenges.plonk_zeta),
        &proof.openings.to_fri_openings(),
        &challenges.fri_challenges,
        merkle_caps,
        &proof.opening_proof,
        &common_data.fri_params,
    )?;

    Ok(())
}
```

`eval_vanishing_poly` 与 `prove` 中计算类似，但只在一个点上计算:

```rust
/// Evaluate the vanishing polynomial at `x`. In this context, the vanishing polynomial is a random
/// linear combination of gate constraints, plus some other terms relating to the permutation
/// argument. All such terms should vanish on `H`.
pub(crate) fn eval_vanishing_poly<F: RichField + Extendable<D>, const D: usize>(
    common_data: &CommonCircuitData<F, D>,
    x: F::Extension,
    vars: EvaluationVars<F, D>,
    local_zs: &[F::Extension],
    next_zs: &[F::Extension],
    local_lookup_zs: &[F::Extension],
    next_lookup_zs: &[F::Extension],
    partial_products: &[F::Extension],
    s_sigmas: &[F::Extension],
    betas: &[F],
    gammas: &[F],
    alphas: &[F],
    deltas: &[F],
) -> Vec<F::Extension> {
    let has_lookup = common_data.num_lookup_polys != 0;
    let max_degree = common_data.quotient_degree_factor;
    let num_prods = common_data.num_partial_products;

    let constraint_terms = evaluate_gate_constraints::<F, D>(common_data, vars);

    let lookup_selectors = &vars.local_constants[common_data.selectors_info.num_selectors()
        ..common_data.selectors_info.num_selectors() + common_data.num_lookup_selectors];

    // The L_0(x) (Z(x) - 1) vanishing terms.
    let mut vanishing_z_1_terms = Vec::new();

    // The terms checking the lookup constraints, if any.
    let mut vanishing_all_lookup_terms = if has_lookup {
        let num_sldc_polys = common_data.num_lookup_polys - 1;
        Vec::with_capacity(
            common_data.config.num_challenges * (4 + common_data.luts.len() + 2 * num_sldc_polys),
        )
    } else {
        Vec::new()
    };

    // The terms checking the partial products.
    let mut vanishing_partial_products_terms = Vec::new();

    let l_0_x = plonk_common::eval_l_0(common_data.degree(), x);

    for i in 0..common_data.config.num_challenges {
        let z_x = local_zs[i];
        let z_gx = next_zs[i];
        vanishing_z_1_terms.push(l_0_x * (z_x - F::Extension::ONE));

        if has_lookup {
            let cur_local_lookup_zs = &local_lookup_zs
                [common_data.num_lookup_polys * i..common_data.num_lookup_polys * (i + 1)];
            let cur_next_lookup_zs = &next_lookup_zs
                [common_data.num_lookup_polys * i..common_data.num_lookup_polys * (i + 1)];

            let cur_deltas = &deltas[NUM_COINS_LOOKUP * i..NUM_COINS_LOOKUP * (i + 1)];

            let lookup_constraints = check_lookup_constraints(
                common_data,
                vars,
                cur_local_lookup_zs,
                cur_next_lookup_zs,
                lookup_selectors,
                cur_deltas.try_into().unwrap(),
            );

            vanishing_all_lookup_terms.extend(lookup_constraints);
        }

        let numerator_values = (0..common_data.config.num_routed_wires)
            .map(|j| {
                let wire_value = vars.local_wires[j];
                let k_i = common_data.k_is[j];
                let s_id = x.scalar_mul(k_i);
                wire_value + s_id.scalar_mul(betas[i]) + gammas[i].into()
            })
            .collect::<Vec<_>>();
        let denominator_values = (0..common_data.config.num_routed_wires)
            .map(|j| {
                let wire_value = vars.local_wires[j];
                let s_sigma = s_sigmas[j];
                wire_value + s_sigma.scalar_mul(betas[i]) + gammas[i].into()
            })
            .collect::<Vec<_>>();

        // The partial products considered for this iteration of `i`.
        let current_partial_products = &partial_products[i * num_prods..(i + 1) * num_prods];
        // Check the quotient partial products.
        let partial_product_checks = check_partial_products(
            &numerator_values,
            &denominator_values,
            current_partial_products,
            z_x,
            z_gx,
            max_degree,
        );
        vanishing_partial_products_terms.extend(partial_product_checks);
    }

    let vanishing_terms = [
        vanishing_z_1_terms,
        vanishing_partial_products_terms,
        vanishing_all_lookup_terms,
        constraint_terms,
    ]
    .concat();

    let alphas = &alphas.iter().map(|&a| a.into()).collect::<Vec<_>>();
    plonk_common::reduce_with_powers_multi(&vanishing_terms, alphas)
}
```

#### Fri Proof



```rust
pub fn verify_fri_proof<
    F: RichField + Extendable<D>,
    C: GenericConfig<D, F = F>,
    const D: usize,
>(
    instance: &FriInstanceInfo<F, D>,
    openings: &FriOpenings<F, D>,
    challenges: &FriChallenges<F, D>,
    initial_merkle_caps: &[MerkleCap<F, C::Hasher>],
    proof: &FriProof<F, C::Hasher, D>,
    params: &FriParams,
) -> Result<()> {
    validate_fri_proof_shape::<F, C, D>(proof, instance, params)?;

    // Size of the LDE domain.
    let n = params.lde_size();

    // Check PoW.
    fri_verify_proof_of_work(challenges.fri_pow_response, &params.config)?;

    // Check that parameters are coherent.
    ensure!(
        params.config.num_query_rounds == proof.query_round_proofs.len(),
        "Number of query rounds does not match config."
    );

    let precomputed_reduced_evals =
        PrecomputedReducedOpenings::from_os_and_alpha(openings, challenges.fri_alpha);
    for (&x_index, round_proof) in challenges
        .fri_query_indices
        .iter()
        .zip(&proof.query_round_proofs)
    {
        fri_verifier_query_round::<F, C, D>(
            instance,
            challenges,
            &precomputed_reduced_evals,
            initial_merkle_caps,
            proof,
            x_index,
            n,
            round_proof,
            params,
        )?;
    }

    Ok(())
}
```



```rust
pub(crate) fn validate_fri_proof_shape<F, C, const D: usize>(
    proof: &FriProof<F, C::Hasher, D>,
    instance: &FriInstanceInfo<F, D>,
    params: &FriParams,
) -> anyhow::Result<()>
where
    F: RichField + Extendable<D>,
    C: GenericConfig<D, F = F>,
{
    let FriProof {
        commit_phase_merkle_caps,
        query_round_proofs,
        final_poly,
        pow_witness: _pow_witness,
    } = proof;

    let cap_height = params.config.cap_height;
    for cap in commit_phase_merkle_caps {
        ensure!(cap.height() == cap_height);
    }

    for query_round in query_round_proofs {
        let FriQueryRound {
            initial_trees_proof,
            steps,
        } = query_round;

        ensure!(initial_trees_proof.evals_proofs.len() == instance.oracles.len());
        for ((leaf, merkle_proof), oracle) in initial_trees_proof
            .evals_proofs
            .iter()
            .zip(&instance.oracles)
        {
            ensure!(leaf.len() == oracle.num_polys + salt_size(oracle.blinding && params.hiding));
            ensure!(merkle_proof.len() + cap_height == params.lde_bits());
        }

        ensure!(steps.len() == params.reduction_arity_bits.len());
        let mut codeword_len_bits = params.lde_bits();
        for (step, arity_bits) in steps.iter().zip(&params.reduction_arity_bits) {
            let FriQueryStep {
                evals,
                merkle_proof,
            } = step;

            let arity = 1 << arity_bits;
            codeword_len_bits -= arity_bits;

            ensure!(evals.len() == arity);
            ensure!(merkle_proof.len() + cap_height == codeword_len_bits);
        }
    }

    ensure!(final_poly.len() == params.final_poly_len());

    Ok(())
}
```



```rust
pub(crate) fn fri_verify_proof_of_work<F: RichField + Extendable<D>, const D: usize>(
    fri_pow_response: F,
    config: &FriConfig,
) -> Result<()> {
    ensure!(
        fri_pow_response.to_canonical_u64().leading_zeros()
            >= config.proof_of_work_bits + (64 - F::order().bits()) as u32,
        "Invalid proof of work witness."
    );

    Ok(())
}
```



```rust
pub(crate) fn from_os_and_alpha(openings: &FriOpenings<F, D>, alpha: F::Extension) -> Self {
    let reduced_openings_at_point = openings
        .batches
        .iter()
        .map(|batch| ReducingFactor::new(alpha).reduce(batch.values.iter()))
        .collect();
    Self {
        reduced_openings_at_point,
    }
}
```



```rust
fn fri_verifier_query_round<
    F: RichField + Extendable<D>,
    C: GenericConfig<D, F = F>,
    const D: usize,
>(
    instance: &FriInstanceInfo<F, D>,
    challenges: &FriChallenges<F, D>,
    precomputed_reduced_evals: &PrecomputedReducedOpenings<F, D>,
    initial_merkle_caps: &[MerkleCap<F, C::Hasher>],
    proof: &FriProof<F, C::Hasher, D>,
    mut x_index: usize,
    n: usize,
    round_proof: &FriQueryRound<F, C::Hasher, D>,
    params: &FriParams,
) -> Result<()> {
    fri_verify_initial_proof::<F, C::Hasher>(
        x_index,
        &round_proof.initial_trees_proof,
        initial_merkle_caps,
    )?;
    // `subgroup_x` is `subgroup[x_index]`, i.e., the actual field element in the domain.
    let log_n = log2_strict(n);
    let mut subgroup_x = F::MULTIPLICATIVE_GROUP_GENERATOR
        * F::primitive_root_of_unity(log_n).exp_u64(reverse_bits(x_index, log_n) as u64);

    // old_eval is the last derived evaluation; it will be checked for consistency with its
    // committed "parent" value in the next iteration.
    let mut old_eval = fri_combine_initial::<F, C, D>(
        instance,
        &round_proof.initial_trees_proof,
        challenges.fri_alpha,
        subgroup_x,
        precomputed_reduced_evals,
        params,
    );

    for (i, &arity_bits) in params.reduction_arity_bits.iter().enumerate() {
        let arity = 1 << arity_bits;
        let evals = &round_proof.steps[i].evals;

        // Split x_index into the index of the coset x is in, and the index of x within that coset.
        let coset_index = x_index >> arity_bits;
        let x_index_within_coset = x_index & (arity - 1);

        // Check consistency with our old evaluation from the previous round.
        ensure!(evals[x_index_within_coset] == old_eval);

        // Infer P(y) from {P(x)}_{x^arity=y}.
        old_eval = compute_evaluation(
            subgroup_x,
            x_index_within_coset,
            arity_bits,
            evals,
            challenges.fri_betas[i],
        );

        verify_merkle_proof_to_cap::<F, C::Hasher>(
            flatten(evals),
            coset_index,
            &proof.commit_phase_merkle_caps[i],
            &round_proof.steps[i].merkle_proof,
        )?;

        // Update the point x to x^arity.
        subgroup_x = subgroup_x.exp_power_of_2(arity_bits);

        x_index = coset_index;
    }

    // Final check of FRI. After all the reductions, we check that the final polynomial is equal
    // to the one sent by the prover.
    ensure!(
        proof.final_poly.eval(subgroup_x.into()) == old_eval,
        "Final polynomial evaluation is invalid."
    );

    Ok(())
}
```



```rust
fn fri_verify_initial_proof<F: RichField, H: Hasher<F>>(
    x_index: usize,
    proof: &FriInitialTreeProof<F, H>,
    initial_merkle_caps: &[MerkleCap<F, H>],
) -> Result<()> {
    for ((evals, merkle_proof), cap) in proof.evals_proofs.iter().zip(initial_merkle_caps) {
        verify_merkle_proof_to_cap::<F, H>(evals.clone(), x_index, cap, merkle_proof)?;
    }

    Ok(())
}

/// Verifies that the given leaf data is present at the given index in the Merkle tree with the
/// given cap.
pub fn verify_merkle_proof_to_cap<F: RichField, H: Hasher<F>>(
    leaf_data: Vec<F>,
    leaf_index: usize,
    merkle_cap: &MerkleCap<F, H>,
    proof: &MerkleProof<F, H>,
) -> Result<()> {
    let mut index = leaf_index;
    let mut current_digest = H::hash_or_noop(&leaf_data);
    for &sibling_digest in proof.siblings.iter() {
        let bit = index & 1;
        index >>= 1;
        current_digest = if bit == 1 {
            H::two_to_one(sibling_digest, current_digest)
        } else {
            H::two_to_one(current_digest, sibling_digest)
        }
    }
    ensure!(
        current_digest == merkle_cap.0[index],
        "Invalid Merkle proof."
    );

    Ok(())
}
```



```rust
pub(crate) fn fri_combine_initial<
    F: RichField + Extendable<D>,
    C: GenericConfig<D, F = F>,
    const D: usize,
>(
    instance: &FriInstanceInfo<F, D>,
    proof: &FriInitialTreeProof<F, C::Hasher>,
    alpha: F::Extension,
    subgroup_x: F,
    precomputed_reduced_evals: &PrecomputedReducedOpenings<F, D>,
    params: &FriParams,
) -> F::Extension {
    assert!(D > 1, "Not implemented for D=1.");
    let subgroup_x = F::Extension::from_basefield(subgroup_x);
    let mut alpha = ReducingFactor::new(alpha);
    let mut sum = F::Extension::ZERO;

    for (batch, reduced_openings) in instance
        .batches
        .iter()
        .zip(&precomputed_reduced_evals.reduced_openings_at_point)
    {
        let FriBatchInfo { point, polynomials } = batch;
        let evals = polynomials
            .iter()
            .map(|p| {
                let poly_blinding = instance.oracles[p.oracle_index].blinding;
                let salted = params.hiding && poly_blinding;
                proof.unsalted_eval(p.oracle_index, p.polynomial_index, salted)
            })
            .map(F::Extension::from_basefield);
        let reduced_evals = alpha.reduce(evals);
        let numerator = reduced_evals - *reduced_openings;
        let denominator = subgroup_x - *point;
        sum = alpha.shift(sum);
        sum += numerator / denominator;
    }

    sum
}
```



```rust
/// Computes P'(x^arity) from {P(x*g^i)}_(i=0..arity), where g is a `arity`-th root of unity
/// and P' is the FRI reduced polynomial.
pub(crate) fn compute_evaluation<F: Field + Extendable<D>, const D: usize>(
    x: F,
    x_index_within_coset: usize,
    arity_bits: usize,
    evals: &[F::Extension],
    beta: F::Extension,
) -> F::Extension {
    let arity = 1 << arity_bits;
    debug_assert_eq!(evals.len(), arity);

    let g = F::primitive_root_of_unity(arity_bits);

    // The evaluation vector needs to be reordered first.
    let mut evals = evals.to_vec();
    reverse_index_bits_in_place(&mut evals);
    let rev_x_index_within_coset = reverse_bits(x_index_within_coset, arity_bits);
    let coset_start = x * g.exp_u64((arity - rev_x_index_within_coset) as u64);
    // The answer is gotten by interpolating {(x*g^i, P(x*g^i))} and evaluating at beta.
    let points = g
        .powers()
        .map(|y| (coset_start * y).into())
        .zip(evals)
        .collect::<Vec<_>>();
    let barycentric_weights = barycentric_weights(&points);
    interpolate(&points, beta, &barycentric_weights)
}
```





### Others



#### Poseidon Hash





#### Merkle Tree



```rust
pub fn index(&self, num_wires: usize, degree: usize) -> usize {
    match self {
        Target::Wire(Wire { row, column }) => row * num_wires + column,
        Target::VirtualTarget { index } => degree * num_wires + index,
    }
}
```



