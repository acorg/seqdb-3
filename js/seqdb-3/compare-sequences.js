
function main()
{
    console.log("compare_sequences_data", compare_sequences_data);
    show_most_frequent_per_group(document.querySelector("#compare-sequences .most-frequent-per-group"))
    show_frequency_per_group(document.querySelector("#compare-sequences .frequency-per-group"))
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

    const add_seqence = function(tr, seq_s, master) {
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
        tr.innerHTML = `<td colspan="2"></td>`;
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
                td.innerHTML = ".";
            }
            else
                td.innerHTML = ".";
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
    for (let group of compare_sequences_data.groups) {
        tab1.appendChild(add_ruler());
        group.seq.forEach(function(id_seq, index) {
            const tr = document.createElement("tr");
            if (index == 0)
                tr.innerHTML = `<td class="group-name" rowspan="${group.seq.length}">${group.name}</td><td class="seq-id">${id_seq.id}</td>`;
            else
                tr.innerHTML = `<td class="seq-id">${id_seq.id}</td>`;
            add_seqence(tr, id_seq.seq, id_seq.id == compare_sequences_data.groups[0].seq[0].id ? "" : compare_sequences_data.groups[0].seq[0].seq);
            tab1.appendChild(tr);
        });
    }
    tab1.appendChild(add_ruler());
    div.appendChild(tab1);
}

// --------------------------------------------------------------------------------

function show_most_frequent_per_group(div)
{
    const title = document.createElement("p");
    title.classList.add("title");
    title.innerHTML = "Most frequent per group";
    div.appendChild(title);


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
