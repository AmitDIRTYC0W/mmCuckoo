const std = @import("std");
const testing = std.testing;

pub const ParameterDefinition = struct {
    minimum: f16,
    range: f16,
};

fn multidimLevySample(
    comptime parameters: []ParameterDefinition,
    step: f16,
    rand: std.rand.Random,
) @Vector(parameters.len, f16) {
    // Generate a random vector
    var random_arr: [parameters.len]f16 = undefined;

    for (&random_arr) |*v| {
        v.* = @floatCast(f16, rand.floatNorm(f32));
    }
    const random_vector: @Vector(parameters.len, f16) = random_arr;

    // Normalise the random vector
    const magnitude = @sqrt(@reduce(.Add, random_vector * random_vector));

    // Choose a target magitude according to the normal distribution (magnitude = c / n^2)
    const n = @floatCast(f16, rand.floatNorm(f32));

    // This will transform the target magnitude to the L`evy distribution. We multiply by a factor
    // of 3 so that the mode of the step size's distribution will be equal to step.
    const scale = step * 3 / (@sqrt(n) * magnitude);

    comptime var ranges: @Vector(parameters.len, f16) = init: {
        var r: @Vector(parameters.len, f16) = undefined;

        for (parameters, 0..) |parameter, i| {
            r[i] = parameter.range;
        }

        break :init r;
    };

    return random_vector * @splat(parameters.len, scale) * ranges;
}

fn Solution(
    comptime parameters: []ParameterDefinition,
    comptime evaluate: fn (parameters: @Vector(parameters.len, f16)) f16,
) type {
    return struct {
        parameters: @Vector(parameters.len, f16),
        quality: f16,

        const Self = @This();

        fn generate(rand_: std.rand.Random) Self {
            var values: [parameters.len]f16 = undefined;
            for (&values, 0..) |*value, i| {
                value.* = parameters[i].minimum + parameters[i].range * @floatCast(f16, rand_.float(f32));
            }

            return .{
                .parameters = values,
                .quality = evaluate(values),
            };
        }

        fn randomStep(self: Self, size: f16, rand_: std.rand.Random) @Vector(parameters.len, f16) {
            while (true) {
                const mod_parameters = self.parameters + multidimLevySample(parameters, size, rand_);

                // Ensure the new parameters are within boundaries
                comptime var minimums: @Vector(parameters.len, f16) = undefined;
                comptime var maximums: @Vector(parameters.len, f16) = undefined;
                comptime {
                    for (parameters, 0..) |parameter, i| {
                        minimums[i] = parameter.minimum;
                        maximums[i] = parameter.minimum + parameter.range;
                    }
                }

                const greater = mod_parameters >= minimums;
                const lower = mod_parameters <= maximums;
                if (@reduce(.And, greater) and @reduce(.And, lower)) {
                    return mod_parameters;
                }
            }
        }

        fn worseThan(_: void, lhs: Self, rhs: Self) bool {
            return lhs.quality < rhs.quality;
        }

        pub fn search(
            population: u16,
            abandon: u16,
            generations: u16,
            step: f16,
            rand_: std.rand.Random,
            allocator: std.mem.Allocator,
        ) !Self {
            if (population < abandon) {
                return error.InvalidParameters;
            }

            // Automatically determine the proper step size
            const step_size = step / @sqrt(@intToFloat(f16, parameters.len * generations));

            // Generate the initial population
            var solutions = try allocator.alloc(Self, population);
            defer allocator.free(solutions);
            for (solutions) |*solution| {
                solution.* = Self.generate(rand_);
            }

            for (0..generations) |_| {
                // Pick a cuckoo randomly
                const i = rand_.uintLessThan(u16, population);
                const j = rand_.uintLessThan(u16, population);

                // Generate a solution by L´evy flights
                const mod_parameters = solutions[i].randomStep(step_size, rand_);

                // Evaluate it
                const mod_quality = evaluate(mod_parameters);

                if (mod_quality > solutions[j].quality) {
                    // Replace solution j with the modified solution
                    solutions[j].parameters = mod_parameters;
                    solutions[j].quality = mod_quality;

                    // Abandon the worst solutions
                    std.sort.heap(Self, solutions, {}, Self.worseThan);

                    for (solutions[0..abandon]) |*solution| {
                        solution.* = Self.generate(rand_);
                    }

                    {
                        const best = solutions[solutions.len - 1];
                        std.debug.print("Best solution: {} of quality {}\n", .{ best.parameters, best.quality });
                    }
                }
            }

            // Return the best solution
            return solutions[solutions.len - 1];
        }
    };
}

fn himmelblau(xy: @Vector(2, f16)) f16 {
    const a = xy[0] * xy[0] + xy[1] - 11;
    const b = xy[0] + xy[1] * xy[1] - 7;
    return -(a * a + b * b);
}

test "Cuckoo search of Himmelblau's function" {
    comptime var parameters = init: {
        var arr: [2]ParameterDefinition = undefined;
        arr[0] = .{ .minimum = -6, .range = 12 };
        arr[1] = .{ .minimum = -6, .range = 12 };
        break :init arr;
    };

    var prng = std.rand.DefaultPrng.init(5);

    const HimmelblauSolution = Solution(&parameters, himmelblau);
    const bestSolution = try HimmelblauSolution.search(
        50,
        20,
        300,
        0.03,
        prng.random(),
        std.testing.allocator,
    );

    try testing.expectApproxEqAbs(bestSolution.quality, 0, 0.01);
}
