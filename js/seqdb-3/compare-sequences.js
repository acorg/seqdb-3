
function main()
{
    console.log("compare_sequences_data", compare_sequences_data);
    show_most_frequent_per_group(document.querySelector("#compare-sequences .most-frequent-per-group"));
    show_frequency_per_group(document.querySelector("#compare-sequences .frequency-per-group"));
    show_positions_with_diversity(document.querySelector("#compare-sequences .positions-with-diversity"));
    show_full_sequences(document.querySelector("#compare-sequences .full-sequences"));
}

// --------------------------------------------------------------------------------

function show_full_sequences(div)
{
    const find_master = function(group) {
        let index = 0;
        for (const seq of group.seq) {
            if (Object.keys(group.pos1).every(function(pos1) { return group.pos1[pos1][0].a == seq.seq[pos1 - 1]; }))
                return index;
            ++index;
        }
        console.warn("master not found");
        return 0;
    };

    // move master sequence (the one having most frequent aas at all positions) to the first element of the group
    const rearrange_group = function(group) {
        const master_index = find_master(group);
        if (master_index > 0)
            group.seq.splice(0, 0, group.seq.splice(master_index, 1)[0]);
    };

    const add_sequence = function(tr, seq_s, master) {
        const seq = [...seq_s];
        seq.forEach(function(aa, pos0) {
            const aa_td = document.createElement("td");
            aa_td.classList.add(`aa${aa}`);
            aa_td.classList.add("aa");
            if (pos0 % 10 == 9)
                aa_td.classList.add("sep-left-zero");
            else if (pos0 % 10 == 4)
                aa_td.classList.add("sep-left-six");
            if (master.length > pos0 && aa == master[pos0])
                aa_td.innerHTML = '&#xB7;';
            else
                aa_td.innerHTML = aa;
            tr.appendChild(aa_td);
        });
    };

    const add_ruler = function() {
        const tr = document.createElement("tr");
        tr.classList.add("aa-ruler");
        tr.innerHTML = `<td colspan="2"></td>`; // group-name + sequence-name
        for (let pos1 = 1; pos1 < 550; ++pos1) {
            const td = document.createElement("td");
            if (pos1 % 10 == 0) {
                td.classList.add("sep-left-zero");
                let pos_s = "" + pos1;
                td.innerHTML = pos_s;
                td.setAttribute("colspan", pos_s.length);
                pos1 += pos_s.length - 1;
            }
            else if (pos1 % 10 == 5) {
                td.classList.add("sep-left-six");
                //td.innerHTML = ".";
            }
            // else
            //     td.innerHTML = ".";
            tr.appendChild(td);
        }
        return tr;
    };

    // ----------------------------------------------------------------------

    const title = document.createElement("p");
    title.classList.add("title");
    title.innerHTML = "Full sequences";
    div.appendChild(title);

    rearrange_group(compare_sequences_data.groups[0]);

    const tab1 = document.createElement("table");
    tab1.appendChild(add_ruler());
    for (let group of compare_sequences_data.groups) {
        group.seq.forEach(function(id_seq, index) {
            const is_master = id_seq.id == compare_sequences_data.groups[0].seq[0].id;
            if (index == 0 && !is_master) {
                const tr_space = document.createElement("tr");
                tr_space.classList.add("group-space");
                tr_space.innerHTML = `<td colspan="${id_seq.seq.length + 2}"></td>`;
                tab1.appendChild(tr_space);
            }
            const tr = document.createElement("tr");
            if (index == 0)
                tr.innerHTML = `<td class="group-name" rowspan="${group.seq.length}">${group.name}</td><td class="seq-id">${id_seq.id}</td>`;
            else
                tr.innerHTML = `<td class="seq-id">${id_seq.id}</td>`;
            add_sequence(tr, id_seq.seq, is_master ? "" : compare_sequences_data.groups[0].seq[0].seq);
            tab1.appendChild(tr);
        });
    }
    tab1.appendChild(add_ruler());
    div.appendChild(tab1);
}

// --------------------------------------------------------------------------------

function position_ruler(positions, initial_tds)
{
        const tr = document.createElement("tr");
        tr.classList.add("position-ruler");
        tr.innerHTML = `<td colspan="${initial_tds}"></td>`; // group-name + sequence-name
        positions.forEach(function(pos1, index) {
            const td = document.createElement("td");
            if (index > 0 && (index % 10 == 0 || index % 10 == 5))
                td.classList.add("sep-left-six");
            td.innerHTML = "" + pos1;
            tr.appendChild(td);
        });
        return tr;
}

// --------------------------------------------------------------------------------

function show_positions_with_diversity(div)
{
    const add_ruler = function() { return position_ruler(compare_sequences_data.pos1, 2); };

    const add_sequence = function(tr, seq, master) {
        compare_sequences_data.pos1.forEach(function(pos1, index) {
            const pos0 = pos1 - 1;
            const aa = seq[pos0];
            const aa_td = document.createElement("td");
            aa_td.classList.add(`aa${aa}`);
            aa_td.classList.add("aa");
            if (index > 0 && (index % 10 == 0 || index % 10 == 5))
                aa_td.classList.add("sep-left-six");
            if (master.length > pos0 && aa == master[pos0])
                aa_td.innerHTML = '&#xB7;';
            else
                aa_td.innerHTML = aa;
            tr.appendChild(aa_td);
        });
    };

    // ----------------------------------------------------------------------

    const title = document.createElement("p");
    title.classList.add("title");
    title.innerHTML = "Positions with diversity";
    div.appendChild(title);

    const tab1 = document.createElement("table");
    tab1.appendChild(add_ruler());
    for (let group of compare_sequences_data.groups) {
        group.seq.forEach(function(id_seq, index) {
            const is_master = id_seq.id == compare_sequences_data.groups[0].seq[0].id;
            if (index == 0 && !is_master) {
                const tr_space = document.createElement("tr");
                tr_space.classList.add("group-space");
                tr_space.innerHTML = `<td colspan="${compare_sequences_data.pos1.length + 2}"></td>`;
                tab1.appendChild(tr_space);
            }
            const tr = document.createElement("tr");
            if (index == 0)
                tr.innerHTML = `<td class="group-name" rowspan="${group.seq.length}">${group.name}</td><td class="seq-id">${id_seq.id}</td>`;
            else
                tr.innerHTML = `<td class="seq-id">${id_seq.id}</td>`;
            add_sequence(tr, id_seq.seq, is_master ? "" : compare_sequences_data.groups[0].seq[0].seq);
            tab1.appendChild(tr);
        });
    }
    tab1.appendChild(add_ruler());
    div.appendChild(tab1);
}

// --------------------------------------------------------------------------------

function show_most_frequent_per_group(div)
{
    const add_ruler = function() { return position_ruler(compare_sequences_data.pos1, 2); }

    // ----------------------------------------------------------------------
    
    const title = document.createElement("p");
    title.classList.add("title");
    title.innerHTML = "Most frequent per group";
    div.appendChild(title);

    const tab1 = document.createElement("table");
    tab1.appendChild(add_ruler());
    for (let group of compare_sequences_data.groups) {
    }
    div.appendChild(tab1);

}

// --------------------------------------------------------------------------------

function show_frequency_per_group(div)
{
    const title = document.createElement("p");
    title.classList.add("title");
    title.innerHTML = "Frequency per group";
    div.appendChild(title);
}

// --------------------------------------------------------------------------------

document.addEventListener("DOMContentLoaded", main);
